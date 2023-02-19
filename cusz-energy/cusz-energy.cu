/**
 * @file capi.cu
 * @author Jiannan Tian
 * @brief
 * @version 0.3
 * @date 2022-05-06
 *
 * (C) 2022 by Washington State University, Argonne National Laboratory
 *
 */

#include <thread>
#include <cstdlib>
#include <string>
#include "../cusz-latest/include/cusz.h"
#include "../cusz-latest/include/cuszapi.hh"

#include "cli/quality_viewer.hh"
#include "cli/timerecord_viewer.hh"
#include <assert.h>
extern "C" {
    #include "papi.h"
}

template <typename T>
void f(std::string fname)
{
    int dev {};
    cudaGetDevice(&dev);
    cudaSetDevice(dev);
    std::string const filename = {"stats.csv"};
    nvmlClass nvml(dev, filename);
    
    /* For demo, we use 3600x1800 CESM data. */
    auto len = 3600 * 1800;

    cusz_header header;
    uint8_t*    exposed_compressed;
    uint8_t*    compressed;
    size_t      compressed_len;

    T *d_uncompressed, *h_uncompressed;
    T *d_decompressed, *h_decompressed;

    /* cuSZ requires a 3% overhead on device (not required on host). */
    size_t uncompressed_memlen = len * 1.03;
    size_t decompressed_memlen = uncompressed_memlen;

    /* code snippet for looking at the device array easily */
    auto peek_devdata = [](T* d_arr, size_t num = 20) {
        thrust::for_each(thrust::device, d_arr, d_arr + num, [=] __device__ __host__(const T i) { printf("%f\t", i); });
        printf("\n");
    };

    // clang-format off
    cudaMalloc(     &d_uncompressed, sizeof(T) * uncompressed_memlen );
    cudaMallocHost( &h_uncompressed, sizeof(T) * len );
    cudaMalloc(     &d_decompressed, sizeof(T) * decompressed_memlen );
    cudaMallocHost( &h_decompressed, sizeof(T) * len );
    // clang-format on

    /* User handles loading from filesystem & transferring to device. */
    io::read_binary_to_array(fname, h_uncompressed, len);
    cudaMemcpy(d_uncompressed, h_uncompressed, sizeof(T) * len, cudaMemcpyHostToDevice);

    /* a casual peek */
    printf("peeking uncompressed data, 20 elements\n");
    peek_devdata(d_uncompressed, 20);

    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // using default
    // cusz_framework* framework = cusz_default_framework();
    // alternatively
    cusz_framework fw = cusz_framework{
        .pipeline     = Auto,
        .predictor    = cusz_custom_predictor{.type = LorenzoI},
        .quantization = cusz_custom_quantization{.radius = 512},
        .codec        = cusz_custom_codec{.type = Huffman}};
    cusz_framework* framework = &fw;

    // Brace initializing a struct pointer is not supported by all host compilers
    // when nvcc forwards.
    // cusz_framework* framework = new cusz_framework{
    //     .pipeline     = Auto,
    //     .predictor    = cusz_custom_predictor{.type = LorenzoI},
    //     .quantization = cusz_custom_quantization{.radius = 512},
    //     .codec        = cusz_custom_codec{.type = Huffman}};

    cusz_compressor* comp       = cusz_create(framework, FP32);
    cusz_config*     config     = new cusz_config{.eb = 1e-1, .mode = Rel};
    cusz_len         uncomp_len = cusz_len{3600, 1800, 1, 1, 1.03};
    cusz_len         decomp_len = uncomp_len;

    std::thread threadStart(&nvmlClass::getStats, &nvml);

    cusz::TimeRecord compress_timerecord;
    cusz::TimeRecord decompress_timerecord;

    {
        cusz_compress(
            comp, config, d_uncompressed, uncomp_len, &exposed_compressed, &compressed_len, &header,
            (void*)&compress_timerecord, stream);

        /* User can interpret the collected time information in other ways. */
        cusz::TimeRecordViewer::view_compression(&compress_timerecord, len * sizeof(T), compressed_len);

        /* verify header */
        printf("header.%-*s : %x\n", 12, "(addr)", &header);
        printf("header.%-*s : %lu, %lu, %lu\n", 12, "{x,y,z}", header.x, header.y, header.z);
        printf("header.%-*s : %lu\n", 12, "filesize", ConfigHelper::get_filesize(&header));
    }

    std::thread threadKill(&nvmlClass::killThread, &nvml);
    threadStart.join();
    threadKill.join();

    /* If needed, User should perform a memcopy to transfer `exposed_compressed` before `compressor` is destroyed. */
    cudaMalloc(&compressed, compressed_len);
    cudaMemcpy(compressed, exposed_compressed, compressed_len, cudaMemcpyDeviceToDevice);

    {
        cusz_decompress(
            comp, &header, exposed_compressed, compressed_len, d_decompressed, decomp_len,
            (void*)&decompress_timerecord, stream);

        cusz::TimeRecordViewer::view_decompression(&decompress_timerecord, len * sizeof(T));
    }

    /* demo: offline checking (de)compression quality. */
    /* load data again    */ cudaMemcpy(d_uncompressed, h_uncompressed, sizeof(T) * len, cudaMemcpyHostToDevice);
    /* perform evaluation */ cusz::QualityViewer::echo_metric_gpu(d_decompressed, d_uncompressed, len, compressed_len);

    cusz_release(comp);

    cudaFree(compressed);
    // delete compressor;

    cudaStreamDestroy(stream);
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        printf("PROG /path/to/cesm-3600x1800\n");
        exit(0);
    }

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    int EventSet = PAPI_NULL;
    long long *values, *values2;
    int num_events=0;
    int code;
    char event_names[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
    char event_descrs[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
    char units[MAX_powercap_EVENTS][PAPI_MIN_STR_LEN];
    int data_type[MAX_powercap_EVENTS];
    int r,i;
    int retval = 0;
    const PAPI_component_info_t *cmpinfo = NULL;
    PAPI_event_info_t evinfo;
    long long before_time,after_time;
    double elapsed_time;

    assert(PAPI_library_init( PAPI_VER_CURRENT ) == PAPI_VER_CURRENT);
    
    assert(PAPI_create_eventset(&EventSet) == PAPI_OK);

    assert(PAPI_add_named_event( EventSet, "nvml:::NVIDIA_GeForce_GTX_1070_Ti:device_0:power" ) == PAPI_OK);
    values= (long long *)calloc(1,sizeof( long long ) );

    assert(PAPI_start(EventSet) == PAPI_OK);

    cudaEventRecord(start);
    f<float>(std::string(argv[1]));
    cudaEventRecord(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    assert(PAPI_stop(EventSet, values) == PAPI_OK);

    power_W = values[0] / 1e3;
    energy_J = power_W * milliseconds * 1e3;

    printf("POWER = %lfW, ENERGY = %lfJ, TIME = %lfms", power_W, energy_J, milliseconds);
    return 0;
}
