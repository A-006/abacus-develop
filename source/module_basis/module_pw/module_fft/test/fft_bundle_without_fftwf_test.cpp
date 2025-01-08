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
class FftBundleTest : public ::testing::Test
{
    protected:
        FFT_Bundle fft_bundle;
        std::string output;
        void SetUp() override
        {
        }
        void TearDown() override
        {
        }
    };

    TEST_F(FftBundleTest,setfft)
    {
        fft_bundle.setfft("cpu","single");
        EXPECT_EQ(fft_bundle.device,"cpu");
        EXPECT_EQ(fft_bundle.precision,"single");

        fft_bundle.setfft("gpu","double");
        EXPECT_EQ(fft_bundle.device,"gpu");
        EXPECT_EQ(fft_bundle.precision,"double");
    }

    TEST_F(FftBundleTest,initfft)
    {
        fft_bundle.setfft("cpu","single");
        testing::internal::CaptureStdout();
        EXPECT_EXIT(fft_bundle.initfft(16,16,16,7,8,256,16,1,false,true),
                            ::testing::ExitedWithCode(1),"");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output,testing::HasSubstr("complie"));

        fft_bundle.setfft("cpu","double");
        fft_bundle.initfft(16,16,16,7,8,256,16,1,false,true);
        EXPECT_EQ(fft_bundle.float_flag,false);
        EXPECT_EQ(fft_bundle.double_flag,true);

        fft_bundle.setfft("gpu","double");
        fft_bundle.initfft(16,16,16,7,8,256,16,1,false,true);
        EXPECT_EQ(fft_bundle.float_flag,false);
        EXPECT_EQ(fft_bundle.double_flag,true);

        fft_bundle.setfft("gpu","single");
        testing::internal::CaptureStdout();
        EXPECT_EXIT(fft_bundle.initfft(16,16,16,7,8,256,16,1,false,true),
                            ::testing::ExitedWithCode(1),"");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output,testing::HasSubstr("complie"));

    }
}