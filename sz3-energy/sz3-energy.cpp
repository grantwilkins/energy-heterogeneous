/*
sz3-energy.c

Test file to run energy testing on generic sz3
compression.
*/

#include "SZ3/api/sz.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

extern "C" {
	#include "papi.h"
}

#define MAX_powercap_EVENTS 64


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


	/* SZ3 COMPRESSION SETUP*/ 
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    float err_bound = 1e-1;
    char oriFilePath[640], outputFilePath[650], decompFilePath[670];
    char *cfgFile;
    size_t cmpSize, oriSize;
    
    if(argc < 3)
    {
		printf("Test case: testfloat_compress [config_file] [srcFilePath] [dimension sizes...] [error bound]\n");
		printf("Example: testfloat_compress sz.config testfloat_8_8_128.dat 8 8 128 1e-3\n");
		exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    if(argc>=4)
		r1 = atoi(argv[3]);
    if(argc>=5)
		r2 = atoi(argv[4]);
    if(argc>=6)
		r3 = atoi(argv[5]);
    
    err_bound = atof(argv[argc - 1]);

    std::vector<size_t> dim_vec = {r1, r2, r3};
   
    printf("cfgFile=%s\n", cfgFile); 
    
    sprintf(outputFilePath, "%s.sz", oriFilePath);
    sprintf(decompFilePath, "%s.sz.out", oriFilePath);

    SZ::Config conf;
    conf.loadcfg(cfgFile);
    conf.errorBoundMode = SZ::EB_REL; // refer to def.hpp for all supported error bound mode
    conf.relErrorBound = err_bound; // relative error bound defaults to 1e-3
    conf.setDims(dim_vec.begin(), dim_vec.end());

    float *data = new float[conf.num];
    SZ::readfile<float>(oriFilePath, conf.num, data);
    size_t outSize;

    assert(PAPI_start(EventSet_comp) == PAPI_OK);
    SZ::Timer timer(true);
    char *bytes = SZ_compress<float>(conf, data, outSize);
    double compress_time = timer.stop();
    assert(PAPI_stop(EventSet_comp, values_comp) == PAPI_OK);
    
    SZ::writefile(outputFilePath, bytes, outSize);

    printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(float) / outSize);
    printf("compression time = %f\n", compress_time);
    printf("compressed data file = %s\n", outputFilePath);


    assert(PAPI_start(EventSet_comp) == PAPI_OK);
    SZ::Timer timer_decomp(true);
    float *decData = SZ_decompress<float>(conf, bytes, cmpSize);
    double decompress_time = timer_decomp.stop();
    assert(PAPI_stop(EventSet_comp, values_comp) == PAPI_OK);


    SZ::writefile<float>(decompFilePath, decData, conf.num);

     //compute the distortion / compression errors...
    SZ::verify<float>(data, decData, conf.num);
    
    printf("decompression time = %f seconds.\n", decompress_time);
    printf("decompressed file = %s\n", outputFilePath);

    delete[]data;
    delete[]bytes;
    delete[] decData;

    /* COMPRESSION & ENERGY MEASUREMENT */


    // STATS REGION //
    // We print the stats out in CSV format in the following order:
    // 1. Data file
    // 2. Error Bound
    // 3. Compression Ratio
    // 4. Compression Runtime
    // 5. Compression CPU Energy
    // +  Compression DRAM Energy
    // 6. Compression CPU Power
    // +  Compression DRAM Power
    // 7. Compression Throughput
    // 8. Decompression Runtime
    // 9. Decompression CPU Energy
    // +  Decompression DRAM Energy
    // 10. Decompression Power
    // +  Decompression DRAM Power
    // 11. Decompression Throughput



    auto cmpData = SZ::readfile<float>(outputFilePath, cmpSize);
    auto oriData = SZ::readfile<float>(oriFilePath, oriSize);

    double runtime_comp, runtime_decomp; // runtime values

    double energy_cpu_comp = 0.0, energy_cpu_decomp = 0.0,
    	energy_dram_comp = 0.0, energy_dram_decomp = 0.0;

    double comp_ratio = double(oriSize)*sizeof(float)/double(cmpSize);

    double throughput_comp = (double(oriSize)/runtime_comp)/1.0e9; // GB/s
    double throughput_decomp = (double(cmpSize)/runtime_decomp)/1.0e9; // GB/s

    runtime_comp = double(end_comp - start_comp) / double(CLOCKS_PER_SEC); // Runtime in seconds
    energy_cpu_comp = double(values_comp[0] + values_comp[1])/(1.0e6); // Energy of CPU in Joules
    energy_dram_comp = double(values_comp[2] + values_comp[3])/(1.0e6); // Energy of CPU in Joules

    runtime_decomp = double(end_comp - start_comp) / double(CLOCKS_PER_SEC);
	energy_cpu_decomp = double(values_decomp[0] + values_decomp[1])/(1.0e6); // Energy of CPU in Joules
    energy_dram_decomp = double(values_decomp[2] + values_decomp[3])/(1.0e6); // Energy of CPU in Joules

    
    FILE * fp;
   	fp = fopen("data.csv", "a");
    fprintf(fp,"%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
    	oriFilePath,
    	err_bound,
    	comp_ratio,
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
    return 0;
}