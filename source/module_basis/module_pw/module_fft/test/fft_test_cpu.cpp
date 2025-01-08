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
template<typename FFT_BASE, typename... Args>
std::unique_ptr<FFT_BASE> make_unique(Args &&... args)
{
    return std::unique_ptr<FFT_BASE>(new FFT_BASE(std::forward<Args>(args)...));
}

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
        fft_cpu_double->initfft(16,16,16,8,8,256,16,1,false,true);
        fft_cpu_double->setupFFT();

        input_float.resize(4096);
        value =0;
        std::generate(input_float.begin(),input_float.end(),
                     [&value](){return std::complex<float>(value++,value);});

        fft_cpu_float =  make_unique<FFT_CPU<float>>(0);
        fft_cpu_float->initfft(16,16,16,8,8,256,16,1,false,true);
        fft_cpu_float->setupFFT();
    }
    void TearDown() override
    {
        fft_cpu_double->clear();
        fft_cpu_float->clear();
    }
};

TEST_F(FftCpuTest,initfft)
{

    EXPECT_EQ(fft_cpu_double->fftnx,16);
    EXPECT_EQ(fft_cpu_double->fftny,16);
    EXPECT_EQ(fft_cpu_double->fftnxy,256);
    EXPECT_EQ(fft_cpu_double->ns,256);
    EXPECT_EQ(fft_cpu_double->lixy,8);
    EXPECT_EQ(fft_cpu_double->rixy,8);
    EXPECT_EQ(fft_cpu_double->nplane,16);

    EXPECT_EQ(fft_cpu_float->fftnx,16);
    EXPECT_EQ(fft_cpu_float->fftny,16);
    EXPECT_EQ(fft_cpu_float->fftnxy,256);
    EXPECT_EQ(fft_cpu_float->ns,256);
    EXPECT_EQ(fft_cpu_float->lixy,8);
    EXPECT_EQ(fft_cpu_float->rixy,8);
    EXPECT_EQ(fft_cpu_float->nplane,16);
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

    //check float cpu type
    fft_cpu_float->setupFFT();
    EXPECT_EQ(sizeof(fft_cpu_float->c_auxg[0]),8);
    EXPECT_EQ(sizeof(fft_cpu_double->s_rspace),8);

    file = fopen("planf_output.txt", "w");
    fflush(stdout); 
    stdout = file;  

    fftwf_print_plan(fft_cpu_float->planfzbac);
    fftwf_print_plan(fft_cpu_float->planfxbac1);
    fclose(file);
    stdout = fdopen(1, "w"); 

    std::ifstream inputf_file("planf_output.txt");
    getline(inputf_file,output);
    EXPECT_THAT(output,testing::HasSubstr("16"));
}

TEST_F(FftCpuTest,fftzfor)
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

// TEST_F(FftCpuTest,)
class FftBunldeTest : public ::testing::Test {
  protected:
    ModulePW::FFT_Bundle fft_bunlde;
    void SetUp()
    {
    	
    }
    
    void TearDown() {  }
};

TEST_F(FftBunldeTest,setfft)
{
    fft_bunlde.setfft("cpu","single");
    EXPECT_EQ(fft_bunlde.device,"cpu");
    EXPECT_EQ(fft_bunlde.precision,"single");
    
    fft_bunlde.setfft("cpu","double");
    EXPECT_EQ(fft_bunlde.device,"cpu");
    EXPECT_EQ(fft_bunlde.precision,"double");
    
}
}