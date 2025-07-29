#ifndef SFP_CONFIG_H
#define SFP_CONFIG_H

#define SFP_OPENCL
#define SFP_OPENMP
#define SFP_BOOST

#define SFP_TARGET_CL_VERSION 300


#ifdef SFP_OPENCL
#define CL_HPP_TARGET_OPENCL_VERSION SFP_TARGET_CL_VERSION
#endif

namespace sfpConfig
{
    static int pulseParameterPadLength = 15;
}

#endif
