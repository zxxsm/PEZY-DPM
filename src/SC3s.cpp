
/*
 * File: SC3s.cpp
 * Author: Zhuxiaoxong
 * Email: zhuxiaoxiong_cn@163.com
 * Description: PEZY host端调用程序
 */
#ifdef MPI_P
#include <mpi.h>
#endif
#include <cstddef>
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <vector>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "Cdglbvar.hpp"
#include "Cfun.hpp"


// 设置LM空间的大小
// size_t stack_size_per_thread = 2048; // 24 - 8*2 = 8KB(SC3)
size_t stack_size_per_thread = 2048; // 24 - 8   = 16KB(SC3)
size_t GLOBALWORKSIZE = 4096;

inline size_t getFileSize(std::ifstream &file)
{
    file.seekg(0, std::ios::end);
    size_t ret = file.tellg();
    file.seekg(0, std::ios::beg);

    return ret;
}
/*
    cache写入的buffer开关，根据我们的测试，打开后写入性能会有一些提升，对读取影响不大
*/
void SetCacheWriteBuffer(cl::Context& context, size_t idxs, bool enable)
{
  typedef int(*pfnExtSetCacheWriteBuffer)(cl_context, size_t, bool);

  pfnExtSetCacheWriteBuffer clExtSetCacheWriteBuffer =
      (pfnExtSetCacheWriteBuffer)clGetExtensionFunctionAddress("pezy_set_cache_writebuffer");

  if (!clExtSetCacheWriteBuffer) {
    throw cl::Error(-1,
        "clGetExtensionFunctionAddress: Can not get pezy_set_cache_writebuffer");
  }

  if ((clExtSetCacheWriteBuffer(context(), idxs, enable)) != CL_SUCCESS) {
    throw cl::Error(-1, "clExtSetCacheWriteBuffer failed");
  }
}
/*
    设置Local memory（线程栈）大小
*/
void SetStackSizePerThread(cl::Kernel &kernel, size_t size)
{

    // Get stack size modify function.
    typedef CL_API_ENTRY cl_int(CL_API_CALL * pfnPezyExtSetPerThreadStackSize)(cl_kernel kernel, size_t size);
    const auto clExtSetPerThreadStackSize = reinterpret_cast<pfnPezyExtSetPerThreadStackSize>(clGetExtensionFunctionAddress("pezy_set_per_thread_stack_size"));
    if (clExtSetPerThreadStackSize == nullptr)
    {
        throw "pezy_set_per_thread_stack_size not found";
    }

    // Set stack size each thread.
    // Each thread's stack size is 2.5KB(SC2) / 2KB(SC) each thread by default.
    // This function will set stack size 1KB each thread.
    // Each PE has 8 thread.
    // And also PE has 20KB(SC2) / 16KB(SC) Scratch pad.
    // So, we can use
    // 24 - 8 = 16KB(SC3)
    // 20 - 8 = 12KB(SC2)
    // 16 - 8 =  8KB(SC)
    // as a user area.
    clExtSetPerThreadStackSize(kernel(), size);
}
inline void loadFile(std::ifstream &file, std::vector<char> &d, size_t size)
{
    d.resize(size);
    file.read(reinterpret_cast<char *>(d.data()), size);
}
cl::Program createProgram(cl::Context &context, const std::vector<cl::Device> &devices, const std::string &filename)
{
    std::ifstream file;
    file.open(filename, std::ios::in | std::ios::binary);

    if (file.fail())
    {
        throw "can not open kernel file";
    }

    size_t filesize = getFileSize(file);
    std::vector<char> binary_data;
    loadFile(file, binary_data, filesize);

    cl::Program::Binaries binaries;
    binaries.push_back(std::make_pair(&binary_data[0], filesize));

    return cl::Program(context, devices, binaries, nullptr, nullptr);
}

cl::Program createProgram(cl::Context &context, const cl::Device &device, const std::string &filename)
{
    std::vector<cl::Device> devices{device};
    return createProgram(context, devices, filename);
}
void pzc_simulate(Tdpmsim *dpmsim, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc,
                  dTdpmesc *dpmesc, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
                  Tdpmrpt *dpmrpt, Trdpmlph *dpmlph, Trdpmpap *dpmpap,
                  Tdpmmat *dpmmat, Tdpmsub *dpmsub, Trdpmscsr *dpmscsr,
                  Trdpmlam *dpmlam, Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw,
                  Trdpmlbr *dpmlbr, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
                  Tcutoff *cutoff, Txbre *xbre, Tsd *sd0, Trseed *rseed0,
                  Tmpi_para *mpi_para, double *tlaps)
{
    //用于计时
    double time0, time1;

    cl::Context context;
    cl::CommandQueue queue;
    std::vector<cl::Kernel> kernels;
    size_t global_work_size;

    // device_id 需要根据myrank进行计算
    size_t device_id = 0;
    try
    {
        // Get Platform
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        const auto &Platform = platforms[0];

        // Get devices
        std::vector<cl::Device> devices;
        Platform.getDevices(CL_DEVICE_TYPE_DEFAULT, &devices);

        // device_id=1; 
        device_id=mpi_para->myrank;
        if (device_id > devices.size())
        {
            std::cerr << "Invalid device id. Use first device." << std::endl;
            device_id = 0;
        }
        // Use device of device_id.
        const auto &device = devices[device_id];

        std::string device_name;
        device.getInfo(CL_DEVICE_NAME, &device_name);

        size_t global_work_size_[3] = {0};
        device.getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &global_work_size_);

        // Get global work size.
        // sc1-64:  8192 (1024 PEs * 8 threads)
        // sc2   : 15782 (1984 PEs * 8 threads)
        // sc3   : 32768 (4096 PEs * 8 threads)
        // sc3s  :  4096 ( 512 PEs * 8 threads)
        // 可以修改最大工作线程数
        global_work_size = GLOBALWORKSIZE;
        if (device_name.find("PEZY-SC2") != std::string::npos)
        {
            global_work_size = std::min(global_work_size, (size_t)15872);
        }

        std::cout << "Use device : " << device_name << std::endl;
        std::cout << "workitem   : " << global_work_size << std::endl;

        // Create Context.
        context = cl::Context(device);

        // Enable L1 cache write buffer.
        SetCacheWriteBuffer(context, 0, 0);

        // Create CommandQueue.
        queue = cl::CommandQueue(context, device, 0);

        // Create Program.
        // Load compiled binary file and create cl::Program object.
        auto program = createProgram(context, device, "kernel/kernel.pz");

        // Create Kernel.
        // Give kernel name without pzc_ prefix.
        auto kernel = cl::Kernel(program, "dpm");

        // 没有多余的栈空间可以用，后续可以考虑加以利用
        // SetStackSizePerThread(kernel, stack_size_per_thread);


        // create device buffer & write 
        // 这部分公共数据CPU---->SCx用于计算
        auto pz_dpmsim = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmsim));
        auto pz_dpmvox = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmvox));
        auto pz_dpmsrc = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmsrc));
        
        //copy to buffer
        queue.enqueueWriteBuffer(pz_dpmsim, true, 0, sizeof(Tdpmsim), dpmsim);
        queue.enqueueWriteBuffer(pz_dpmvox, true, 0, sizeof(Tdpmvox), dpmvox);
        queue.enqueueWriteBuffer(pz_dpmsrc, true, 0, sizeof(Tdpmsrc), dpmsrc);
        
        auto pz_dpmwck = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmwck));
        auto pz_dpmcmp = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmcmp));

        queue.enqueueWriteBuffer(pz_dpmwck, true, 0, sizeof(Trdpmwck), dpmwck);
        queue.enqueueWriteBuffer(pz_dpmcmp, true, 0, sizeof(Trdpmcmp), dpmcmp);
        
        auto pz_dpmrpt = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmrpt));
        auto pz_dpmlph = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmlph));
        auto pz_dpmpap = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmpap));
        
        queue.enqueueWriteBuffer(pz_dpmrpt, true, 0, sizeof(Tdpmrpt), dpmrpt);
        queue.enqueueWriteBuffer(pz_dpmlph, true, 0, sizeof(Trdpmlph), dpmlph);
        queue.enqueueWriteBuffer(pz_dpmpap, true, 0, sizeof(Trdpmpap), dpmpap);

        auto pz_dpmmat  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmmat));
        auto pz_dpmsub  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tdpmsub));
        auto pz_dpmscsr = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmscsr));

        queue.enqueueWriteBuffer(pz_dpmmat, true, 0, sizeof(Tdpmmat), dpmmat);
        queue.enqueueWriteBuffer(pz_dpmsub, true, 0, sizeof(Tdpmsub), dpmsub);
        queue.enqueueWriteBuffer(pz_dpmscsr, true, 0, sizeof(Trdpmscsr), dpmscsr);

        auto pz_dpmlam  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmlam));
        auto pz_dpmstpw = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmstpw));
        auto pz_dpmscpw = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmscpw));

        queue.enqueueWriteBuffer(pz_dpmlam,  true, 0, sizeof(Trdpmlam ), dpmlam);
        queue.enqueueWriteBuffer(pz_dpmstpw, true, 0, sizeof(Trdpmstpw), dpmstpw);
        queue.enqueueWriteBuffer(pz_dpmscpw, true, 0, sizeof(Trdpmscpw), dpmscpw);

        auto pz_dpmlbr  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trdpmlbr));
        auto pz_xdpmbw  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trxdpmbw));
        auto pz_xdpmq2D = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Txdpmq2D));

        queue.enqueueWriteBuffer(pz_dpmlbr, true, 0, sizeof(Trdpmlbr), dpmlbr);
        queue.enqueueWriteBuffer(pz_xdpmbw, true, 0, sizeof(Trxdpmbw), xdpmbw);
        queue.enqueueWriteBuffer(pz_xdpmq2D, true, 0, sizeof(Txdpmq2D), xdpmq2D);

        auto pz_cutoff  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Tcutoff));
        auto pz_xbre    = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Txbre));
        auto pz_rseed0  = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Trseed));

        queue.enqueueWriteBuffer(pz_cutoff, true, 0, sizeof(Tcutoff), cutoff);
        queue.enqueueWriteBuffer(pz_xbre, true, 0, sizeof(Txbre), xbre);
        queue.enqueueWriteBuffer(pz_rseed0, true, 0, sizeof(Trseed), rseed0);
        
        //需要把SCx和CPU需要的nprtcl计算出来
        //根据不同的CPU，SCx_ratio选择的值也不同，通常由测试来确定
        int32 nprtcl_pezy=_maxhis*SCx_ratio;
        int32 nprtcl_cpu=_maxhis*(1.0-SCx_ratio);

        nprtcl_pezy=Ccomp_np(nprtcl_pezy,mpi_para->myrank,mpi_para->pro_num);
#ifdef DEBUG
        printf("PEZY-SCx,PRO-ID:%d,device_id-%ld,PART-NUM-PEZY:%d\n", mpi_para->myrank,device_id,nprtcl_pezy);
#endif

        //传递Tsd时需要进行种子分割,分割的基本原理为
        //第一个SCx用 [0,global_work_size-1]
        //第二个SCx用 [global_work_size,2*global_work_size-1]
        //……
        //当所有的SCx全部用完了，再考虑CPU的线程使用。
        auto pz_sd0     = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int32)*2*global_work_size);
        queue.enqueueWriteBuffer(pz_sd0, true, 0, sizeof(int32)*2*global_work_size, 
        &(*sd0).Iseed[mpi_para->myrank*2*global_work_size]);
#ifdef DEBUG
        printf("PEZY-SCx,PRO-ID:%d,RAND SEED:%d  %d\n", mpi_para->myrank,
               (*sd0).Iseed[mpi_para->myrank * 2 * global_work_size],
               (*sd0).Iseed[mpi_para->myrank * 2 * global_work_size + 1]);
#endif
        // 注：Tdpmesc是进行规约的参数，CPU程序中可以单个线程/进程单独所有；
        // 但是单个Tdpmesc大小为16.4MB，在SCx中占用内存太大。
        // 这里需要设置开辟多少个空间，且由多个线程共同进行规约
        // 规约后将首个结果传回
        auto pz_dpmesc = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(dTdpmesc));
        
        // Tdpmesc 初始化为0
        cl::Event write_dpmesc;
        queue.enqueueFillBuffer(pz_dpmesc, 0, 0, sizeof(dTdpmesc), nullptr, &write_dpmesc);
        write_dpmesc.wait();
        // 注：Tdpmstck是每个线程私有的参数；SCx中每个线程单独所有；
        // 单个 Tdpmstck 大小为10.6KB , 超过了单个线程的stack 上界；
        // 只能在CPU端对SCx端进行开辟 线程数 个空间；
        auto pz_dpmstck = cl::Buffer(context, CL_MEM_READ_WRITE, global_work_size*sizeof(Trdpmstck));
        
        // 验证参数(可删除)
        auto pz_veryfy = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int32) * 1);
        // Clear dst.
        cl::Event write_verfy;
        queue.enqueueFillBuffer(pz_veryfy, 0, 0, sizeof(int32) * 1, nullptr, &write_verfy);
        //等待写入完成；
        write_verfy.wait();

        //set kernel
        kernel.setArg(0,  pz_dpmsim);
        kernel.setArg(1,  pz_dpmvox);
        kernel.setArg(2,  pz_dpmsrc);
        kernel.setArg(3,  pz_dpmesc);
        kernel.setArg(4,  pz_dpmwck);
        kernel.setArg(5,  pz_dpmcmp);
        kernel.setArg(6,  pz_dpmrpt);
        kernel.setArg(7,  pz_dpmlph);
        kernel.setArg(8,  pz_dpmpap);
        kernel.setArg(9,  pz_dpmmat);
        kernel.setArg(10, pz_dpmsub);
        kernel.setArg(11, pz_dpmscsr);
        kernel.setArg(12, pz_dpmlam);
        kernel.setArg(13, pz_dpmstpw);
        kernel.setArg(14, pz_dpmscpw);
        kernel.setArg(15, pz_dpmlbr);
        kernel.setArg(16, pz_xdpmbw);
        kernel.setArg(17, pz_xdpmq2D);
        kernel.setArg(18, pz_cutoff);
        kernel.setArg(19, pz_xbre);
        kernel.setArg(20, pz_sd0);//已经分了段的种子数
        kernel.setArg(21, pz_rseed0);
        kernel.setArg(22, pz_dpmstck);
        kernel.setArg(23, nprtcl_pezy);//已经计算好的粒子数
        kernel.setArg(24, pz_veryfy);  //用于传回值验证结果

        //对kernel进行记时
#ifdef MPI_P
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        timers(&time0);

        // Run device kernel.
        cl::Event event_dpm;
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(), nullptr, &event_dpm);
        event_dpm.wait();

#ifdef MPI_P
        //数据传回
        dTdpmesc *dpmesc_reduce = (dTdpmesc *)malloc(sizeof(dTdpmesc));
        queue.enqueueReadBuffer(pz_dpmesc, true, 0, sizeof(dTdpmesc), dpmesc_reduce);//规约后将结果传回
        MPI_Reduce(dpmesc_reduce, dpmesc, 2 * nxyz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        queue.enqueueReadBuffer(pz_dpmesc, true, 0, sizeof(dTdpmesc), dpmesc);//规约后将结果传回
#endif

        timers(&time1);   
        *tlaps = time1 - time0;

        // Get dst.
        int32 veryfy;
        queue.enqueueReadBuffer(pz_veryfy, true, 0, sizeof(int32) * 1, &veryfy);
   
        std::cout<<"veryfy="<<veryfy<<std::endl;

        queue.flush();
        queue.finish();
    }
    catch (const cl::Error &e)
    {
        std::stringstream msg;
        msg << "CL Error : " << e.what() << " " << e.err();
        throw std::runtime_error(msg.str());
    }
}
