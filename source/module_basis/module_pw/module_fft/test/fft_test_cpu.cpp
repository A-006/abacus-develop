#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "malloc.h"
#include "fstream"
#define private public
#include "../fft_base.h"
#include "../fft_bundle.h"
#include "../fft_cpu.h"
#include "module_parameter/parameter.h"
#undef private
namespace ModulePW
{

class FftCpuTest : public ::testing::Test
{
protected:
    std::shared_ptr<FFT_CPU<double>> fft_cpu_double;
    std::shared_ptr<FFT_CPU<float>>  fft_cpu_float;
    std::string                      output;
    std::vector<std::complex<double>> input_double;
    std::vector<std::complex<float>>  input_float;
    void SetUp() override
    {
        int value =0;
        input_double.resize(4096);
        std::generate(input_double.begin(),input_double.end(),
                      [&value]{return std::complex<double>(value++,value);});
        fft_cpu_double = make_unique<FFT_CPU<double>>(0);
        fft_cpu_double->initfft(16,16,16,7,8,256,16,1,false,true);
        fft_cpu_double->setupFFT();

    }
    void TearDown() override
    {
        fft_cpu_double->clear();
    }
};

TEST_F(FftCpuTest,initfft)
{

    EXPECT_EQ(fft_cpu_double->fftnx,16);
    EXPECT_EQ(fft_cpu_double->fftny,16);
    EXPECT_EQ(fft_cpu_double->fftnxy,256);
    EXPECT_EQ(fft_cpu_double->ns,256);
    EXPECT_EQ(fft_cpu_double->lixy,7);
    EXPECT_EQ(fft_cpu_double->rixy,8);
    EXPECT_EQ(fft_cpu_double->nplane,16);

}

TEST_F(FftCpuTest,setupFFT)
{
    //cheak double cpu type
    
    EXPECT_EQ(sizeof(fft_cpu_double->z_auxg[0]),16);
    EXPECT_EQ(sizeof(fft_cpu_double->d_rspace),8);

    FILE *file = fopen("plan_output.txt", "w");
    fflush(stdout); 
    stdout = file;  

    fftw_print_plan(fft_cpu_double->planzbac);
    fftw_print_plan(fft_cpu_double->planxbac1);
    fclose(file);
    stdout = fdopen(1, "w"); 

    std::ifstream input_file("plan_output.txt");
    getline(input_file,output);
    input_file.close();
    EXPECT_THAT(output,testing::HasSubstr("16"));

}

TEST_F(FftCpuTest,fftzforAndBac)
{
    fft_cpu_double->fftzfor(input_double.data(),input_double.data());
    fft_cpu_double->fftzbac(input_double.data(),input_double.data());
    
    // the fftw need to normalization
    for (int i=0;i<4096;i++)
    { 
        EXPECT_EQ(input_double[i],
                  std::complex<double>(fft_cpu_double->nplane*i,fft_cpu_double->nplane*i));
    }

}

TEST_F(FftCpuTest,fftxyforAndBac)
{
    fft_cpu_double->fftxyfor(input_double.data(),input_double.data());
    fft_cpu_double->fftxybac(input_double.data(),input_double.data());
    for (int i=0;i<2;i++)
    { 
        EXPECT_EQ(input_double[i],
                  std::complex<double>(fft_cpu_double->fftnxy*i,fft_cpu_double->fftnxy*i));
    }

}

TEST_F(FftCpuTest,clearfft)
{
    //test clean plan
    EXPECT_NE(fft_cpu_double->planzbac,nullptr);
    EXPECT_NE(fft_cpu_double->planzfor,nullptr);
    EXPECT_NE(fft_cpu_double->planxfor1,nullptr);
    EXPECT_NE(fft_cpu_double->planxbac1,nullptr);
    EXPECT_NE(fft_cpu_double->planyfor,nullptr);
    EXPECT_NE(fft_cpu_double->planybac,nullptr);
    fft_cpu_double->cleanFFT();
    EXPECT_EQ(fft_cpu_double->planzbac,nullptr);
    EXPECT_EQ(fft_cpu_double->planzfor,nullptr);
    EXPECT_EQ(fft_cpu_double->planxfor1,nullptr);
    EXPECT_EQ(fft_cpu_double->planxbac1,nullptr);
    EXPECT_EQ(fft_cpu_double->planyfor,nullptr);
    EXPECT_EQ(fft_cpu_double->planybac,nullptr);

}

TEST_F(FftCpuTest,clear)
{
    EXPECT_NE(fft_cpu_double->z_auxg,nullptr);
    EXPECT_NE(fft_cpu_double->z_auxr,nullptr);
    fft_cpu_double->clear();
    EXPECT_EQ(fft_cpu_double->z_auxg,nullptr);
    EXPECT_EQ(fft_cpu_double->z_auxr,nullptr);

    
}
}