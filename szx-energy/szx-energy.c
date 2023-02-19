/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "szx.h"
#include "szx_rw.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "papi.h"
#include <assert.h>

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

#define MAX_powercap_EVENTS 64

void cost_start()
{
    totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    int EventSet_comp = PAPI_NULL, EventSet_decomp = PAPI_NULL;
    long long *values_comp, *values_decomp;
    int num_events=0;
    int code;
    char event_names[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
    char event_descrs[MAX_powercap_EVENTS][PAPI_MAX_STR_LEN];
    char units[MAX_powercap_EVENTS][PAPI_MIN_STR_LEN];
    int data_type[MAX_powercap_EVENTS];
    int r,i;
    int retval = 0;
    clock_t start_comp, end_comp, start_decomp, end_decomp;

    const PAPI_component_info_t *cmpinfo = NULL;
    PAPI_event_info_t evinfo;
    long long before_time,after_time;
    double elapsed_time;

    assert(PAPI_library_init( PAPI_VER_CURRENT ) == PAPI_VER_CURRENT);
    
    assert(PAPI_create_eventset(&EventSet_comp) == PAPI_OK);
    assert(PAPI_create_eventset(&EventSet_decomp) == PAPI_OK);


    assert(PAPI_add_named_event( EventSet_comp, 
        "powercap:::ENERGY_UJ:ZONE0" ) == PAPI_OK);
    /*  
    assert(PAPI_add_named_event( EventSet_comp, 
        "powercap:::ENERGY_UJ:ZONE1" ) == PAPI_OK);
    assert(PAPI_add_named_event( EventSet_comp, 
        "powercap:::ENERGY_UJ:ZONE0_SUBZONE0" ) == PAPI_OK);
    assert(PAPI_add_named_event( EventSet_comp, 
        "powercap:::ENERGY_UJ:ZONE1_SUBZONE0" ) == PAPI_OK);
    */
    assert(PAPI_add_named_event( EventSet_decomp, 
        "powercap:::ENERGY_UJ:ZONE0" ) == PAPI_OK);
    /*
    assert(PAPI_add_named_event( EventSet_decomp, 
        "powercap:::ENERGY_UJ:ZONE1" ) == PAPI_OK);
    assert(PAPI_add_named_event( EventSet_decomp, 
        "powercap:::ENERGY_UJ:ZONE0_SUBZONE0" ) == PAPI_OK);
    assert(PAPI_add_named_event( EventSet_decomp, 
        "powercap:::ENERGY_UJ:ZONE1_SUBZONE0" ) == PAPI_OK);
    */
    values_comp = (long long *)calloc(4, sizeof( long long ));
    values_decomp = (long long *)calloc(4, sizeof( long long ));


    char oriFilePath[640], outputFilePath[645], decompFilePath[670];
    if(argc < 3)
    {
        printf("Usage: testfloat_compress_fastmode3 [srcFilePath] [block size] [err bound]\n");
        printf("Example: testfloat_compress_fastmode3 testfloat_8_8_128.dat 64 1E-3\n");
        exit(0);
    }

    sprintf(oriFilePath, "%s", argv[1]);
    int blockSize = atoi(argv[2]);
    float errBound = atof(argv[3]);

    sprintf(outputFilePath, "%s.szx", oriFilePath);
    sprintf(decompFilePath, "%s.szx.out", oriFilePath);

    int status = 0;
    size_t nbEle;
    size_t oriSize = nbEle*sizeof(float);
    float *ori_data = readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZ_SCES)
    {
        printf("Error: data file %s cannot be read!\n", oriFilePath);
        exit(0);
    }

    size_t outSize;

    assert(PAPI_start(EventSet_comp) == PAPI_OK);
    unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp(ori_data, &outSize, errBound, nbEle, blockSize);
    assert(PAPI_stop(EventSet_comp, values_comp) == PAPI_OK);


    writeByteData(bytes, outSize, outputFilePath, &status);
    if(status != SZ_SCES)
    {
        printf("Error: data file %s cannot be written!\n", outputFilePath);
        exit(0);
    }


    cost_start();
    float *data = NULL;
    //SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&data, nbEle, bytes);
    SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(&data, nbEle, bytes);
    cost_end();

    printf("timecost=%f\n",totalCost);
    writeFloatData_inBytes(data, nbEle, outputFilePath, &status);
    
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be written!\n", outputFilePath);
        exit(0);
    }

    i = 0;
    float Max = 0, Min = 0, diffMax = 0;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = fabs(data[0] - ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < nbEle; i++)
    {
        sum1 += ori_data[i];
        sum2 += data[i];
    }
    double mean1 = sum1/nbEle;
    double mean2 = sum2/nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;

    double maxpw_relerr = 0;
    for (i = 0; i < nbEle; i++)
    {
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];

        float err = fabs(data[i] - ori_data[i]);
    if(ori_data[i]!=0)
    {
        if(fabs(ori_data[i])>1)
            relerr = err/ori_data[i];
        else
            relerr = err;
        if(maxpw_relerr<relerr)
            maxpw_relerr = relerr;
        }

    /*if(err > 1600000)
    {
        printf("i=%zu, ori=%f, dec=%f, diff=%f\n", i, ori_data[i], data[i], err);
        exit(0);
    }*/
    if (diffMax < err)
        diffMax = err;
        prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
        sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
        sum4 += (data[i] - mean2)*(data[i]-mean2);
    sum += err*err;
    }
    double std1 = sqrt(sum3/nbEle);
    double std2 = sqrt(sum4/nbEle);
    double ee = prodSum/nbEle;
    double acEff = ee/std1/std2;

    double mse = sum/nbEle;
    double range = Max - Min;
    double psnr = 20*log10(range)-10*log10(mse);
    double nrmse = sqrt(mse)/range;

    double compressionRatio = 1.0*nbEle*sizeof(float)/outSize;

    printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
    printf ("Max absolute error = %.10f\n", diffMax);
    printf ("Max relative error = %f\n", diffMax/(Max-Min));
    printf ("Max pw relative error = %f\n", maxpw_relerr);
    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
    printf ("acEff=%f\n", acEff);
    printf ("compressionRatio = %f\n", compressionRatio);


    double runtime_comp = totalCost, runtime_decomp = totalCost; // runtime values

    double energy_cpu_comp = 0.0, energy_cpu_decomp = 0.0,
        energy_dram_comp = 0.0, energy_dram_decomp = 0.0;

    //runtime_comp = double(end_comp - start_comp) / double(CLOCKS_PER_SEC); // Runtime in seconds
    energy_cpu_comp = (double)(values_comp[0] + values_comp[1])/(1.0e6); // Energy of CPU in Joules
    energy_dram_comp = (double)(values_comp[2] + values_comp[3])/(1.0e6); // Energy of CPU in Joules

    //runtime_decomp = double(end_comp - start_comp) / double(CLOCKS_PER_SEC);
    energy_cpu_decomp = (double) (values_decomp[0] + values_decomp[1])/(1.0e6); // Energy of CPU in Joules
    energy_dram_decomp = (double)(values_decomp[2] + values_decomp[3])/(1.0e6); // Energy of CPU in Joules

    double throughput_comp = (double) oriSize /(runtime_comp)/1.0e9; // GB/s
    double throughput_decomp = (double) outSize /(runtime_decomp)/1.0e9; // GB/s

    FILE * fp;
    fp = fopen("data.csv", "a");
    fprintf(fp,"%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
        oriFilePath,
        errBound,
        compressionRatio,
        runtime_comp,
        energy_cpu_comp,
        energy_dram_comp,
        energy_cpu_comp/runtime_comp,
        energy_dram_comp/runtime_comp,
        throughput_comp,
        runtime_decomp,
        energy_cpu_decomp,
        energy_dram_decomp,
        energy_cpu_decomp/runtime_decomp,
        energy_dram_decomp/runtime_decomp,
        throughput_decomp);
    fclose(fp);

    free(ori_data);
    free(bytes);
    free(data);

    return 0;
}