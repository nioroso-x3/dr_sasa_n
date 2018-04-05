#include "stdafx.h"
#include "atom_struct.h"
//#define CL_HPP_MINIMUM_OPENCL_VERSION 120
//#define CL_HPP_TARGET_OPENCL_VERSION 120
//#ifdef __APPLE__
//#include <OpenCL/cl.hpp>
//#else
//#include <CL/cl.hpp>
//#endif
//#include <chrono>
#include "SurfaceSolverCL.h"


void
SimpleSolverCL(vector<atom_struct>& pdb,
               vector<float>&      points,
               int                  cl_dev_type){
  uint32 P_SIZE = points.size()/3;
  int err = -999;
  //vector<cl::Platform> all_platforms;
 // cl::Platform::get(&all_platforms);
  //if(all_platforms.size()==0|| cl_dev_type == 0){
    //cerr<<"No platforms found. Using CPU.\n";
    auto* p = points.data();
    uint32 pdbs = pdb.size();
    #pragma omp parallel for schedule(dynamic)
    for(uint32 i = 0; i < pdbs; ++i){
      auto& atom_i = pdb[i];

      if (!atom_i.ACTIVE) continue;

      atom_i.SHELL_BURIED.resize(P_SIZE,false);
      atom_i.ACCESSIBLE_P = P_SIZE;
      auto* b = atom_i.SHELL_BURIED.data();
      auto* C_I = atom_i.COORDS.data();
      auto R_I = atom_i.RADIUS;

      for (uint32 j = 0; j < atom_i.INTERACTION_P.size();++j){
        auto& atom_j = pdb[atom_i.INTERACTION_P[j]];
        if (!atom_j.ACTIVE) continue;
        auto* C_J = atom_j.COORDS.data();
        auto R_J2 = atom_j.RADIUS2;
        uint counter = 0;
        for (uint32 k = 0; k < P_SIZE; ++k){
          if(!b[k]){
            counter++;
            auto k3 = k * 3;
            auto Xj = p[k3]     * R_I + C_I[0] - C_J[0];
            auto Yj = p[k3 + 1] * R_I + C_I[1] - C_J[1];
            auto Zj = p[k3 + 2] * R_I + C_I[2] - C_J[2];
            float dist_j = Xj*Xj + Yj*Yj + Zj*Zj;
            b[k] = (dist_j <= R_J2);
          }
        }
      }
      atom_i.ACCESSIBLE_P -= count(atom_i.SHELL_BURIED.begin(), atom_i.SHELL_BURIED.end(), true);
      atom_i.SASA = atom_i.AREA*((float)atom_i.ACCESSIBLE_P / (float)P_SIZE);
    }
    return;
  //}
 

 /* vector<cl::Device> cpu_devices;
  vector<cl::Device> gpu_devices;
  vector<cl::Device> all_devices;
  for(auto& platform : all_platforms){
    vector<cl::Device> cpud;
    vector<cl::Device> gpud;
    platform.getDevices(CL_DEVICE_TYPE_CPU, &cpud);
    platform.getDevices(CL_DEVICE_TYPE_GPU, &gpud);
    for (auto& dev : cpud){
//       const cl_device_partition_property subDeviceProperties[] = {
//             CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN,
//             CL_DEVICE_AFFINITY_DOMAIN_L3_CACHE,
//             CL_PROPERTIES_LIST_END_EXT,
//             0};
//       vector<cl::Device> subcpu;
//       dev.createSubDevices(subDeviceProperties,&subcpu);
//       for (auto& subdev : subcpu){
//         cpu_devices.push_back(subdev);
//         all_devices.push_back(subdev);
//       }
      cpu_devices.push_back(dev);
      all_devices.push_back(dev);
    }
    for (auto& dev : gpud){
      gpu_devices.push_back(dev);
      all_devices.push_back(dev);
    }
  }
  if (all_devices.size() == 0){
    cerr << "No OpenCL devices found.\n";
    throw -1;
  }
  cl::Device default_device;
  if(cl_dev_type >= 1){
    if(gpu_devices.size()){
      uint32 devN = cl_dev_type - 1;
      if (devN > gpu_devices.size()){
        cerr << "Invalid GPU device\n";
        throw -1;
      }
      else default_device = gpu_devices[devN];
    }  
    else{
      cerr << "No OpenCL GPU devices found, selecting first dev.\n";
      default_device = all_devices[0];
    }
  }
  else{
    if(cpu_devices.size())
      default_device = cpu_devices[0];
    else{
      cerr << "No OpenCL CPU devices found, selecting first dev.\n";
      default_device = all_devices[0];
    }
  }
  
  bool devFP64 = default_device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_fp64") != string::npos;
//   cerr<< "SimpleSolverCL using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<< 
//         (devFP64 ? " FP64\n" : " FP32\n");
  cl::Context context((vector<cl::Device>){default_device});
  cl::Program::Sources sources;
  vector<string> kernels = {
  "#define P_SIZE "+to_string(P_SIZE)+"                                          \n"
  "#pragma OPENCL EXTENSION cl_khr_fp64 : enable                              \n"
  "void kernel simple_solver(global const float*   coord,                    \n"
  "                          global const float*   point,                    \n"
  "                          global const uint*     inter,                    \n"
  "                          global const uint*     interpos,                 \n"
  "                          global uint*           apoint                    \n"
  "){                                                                        \n"
  "  uint gid = get_global_id(0);                                            \n"
  "  uint i_idx = gid / P_SIZE;                                              \n"
  "  uint ci_idx = i_idx*4;                                                  \n"
  "  uint K4 = (gid - i_idx*P_SIZE)*4;                                       \n"
  "  float R_I = coord[ci_idx];                                              \n"
  "  for(uint i = interpos[i_idx]; i < interpos[i_idx + 1]; ++i){            \n"
  "    uint cj_idx = inter[i]*4;                                             \n"
  "    float R_J = coord[cj_idx];                                            \n"
  "    float dX = point[K4]     * R_I + coord[ci_idx + 1] - coord[cj_idx + 1];\n"
  "    float dY = point[K4 + 1] * R_I + coord[ci_idx + 2] - coord[cj_idx + 2];\n"
  "    float dZ = point[K4 + 2] * R_I + coord[ci_idx + 3] - coord[cj_idx + 3];\n"
  "    if((dX*dX + dY*dY + dZ*dZ) <= (R_J*R_J)){                             \n"
  "      atomic_dec(apoint+i_idx);                                           \n"
  "      break;                                                              \n"
  "    }                                                                     \n"
  "  }                                                                       \n"
  "}                                                                         \n",
  "#define P_SIZE "+to_string(P_SIZE)+"                                          \n"
  "void kernel simple_solver(global const float*   coord,                    \n"
  "                          global const float*   point,                    \n"
  "                          global const uint*    inter,                    \n"
  "                          global const uint*    interpos,                 \n"
  "                          global uint*          apoint                    \n"
  "){                                                                        \n"
  "  uint gid = get_global_id(0);                                            \n"
  "  uint i_idx = gid / P_SIZE;                                              \n"
  "  uint ci_idx = i_idx*4;                                                  \n"
  "  uint K4 = (gid - i_idx*P_SIZE)*4;                                       \n"
  "  float R_I = coord[ci_idx];                                              \n"
  "  for(uint i = interpos[i_idx]; i < interpos[i_idx + 1]; ++i){            \n"
  "    uint cj_idx = inter[i]*4;                                             \n"
  "    float R_J = coord[cj_idx];                                            \n"
  "    float dX = point[K4]     * R_I + coord[ci_idx + 1] - coord[cj_idx + 1];\n"
  "    float dY = point[K4 + 1] * R_I + coord[ci_idx + 2] - coord[cj_idx + 2];\n"
  "    float dZ = point[K4 + 2] * R_I + coord[ci_idx + 3] - coord[cj_idx + 3];\n"
  "    if((dX*dX + dY*dY + dZ*dZ) <= (R_J*R_J)){                             \n"
  "      atomic_dec(apoint+i_idx);                                           \n"
  "      break;                                                              \n"
  "    }                                                                     \n"
  "  }                                                                       \n"
  "}                                                                         \n"
  };
  string kernel_code;
  if(devFP64)  kernel_code = kernels[0];
  else         kernel_code = kernels[1];
  sources.push_back({kernel_code.c_str(),kernel_code.length()});
  cl::Program program(context,sources);
  if(program.build({default_device})!=CL_SUCCESS){
    cerr<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
    throw -1;
  }
  //else cerr << "OpenCL kernel built.\n";  
  
  //chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  //cout << "exec\n";
  if(devFP64){
    vector<float> COORD_DATA(pdb.size()*4);
    vector<float> POINT_DATA(P_SIZE*4);
    vector<uint> AREA_COUNT(pdb.size(),P_SIZE);
    vector<uint> INT_DATA;
    vector<uint> INTER_POS;
    //vector<char> AREAP(pdb.size()*P_SIZE,0);
    #pragma omp parallel for
    for (uint32 i = 0; i < pdb.size(); ++i){
      int i4 = i*4;
      COORD_DATA[i4] =     pdb[i].RADIUS;
      COORD_DATA[i4 + 1] = pdb[i].COORDS[0];
      COORD_DATA[i4 + 2] = pdb[i].COORDS[1];
      COORD_DATA[i4 + 3] = pdb[i].COORDS[2];
    }

    INTER_POS.push_back(0);
    for(auto& atom : pdb){
      INTER_POS.push_back(INTER_POS.back()+atom.INTERACTION_P.size());
      //cout << atom.ID << "\t" << atom.STRUCT_TYPE << " | ";
      for(auto& p : atom.INTERACTION_P) {
        //cout << pdb[p].ID<<"|"<<pdb[p].STRUCT_TYPE  << "\t";
        INT_DATA.push_back(p);
      }
      //cout << "\n";
    }
    INT_DATA.shrink_to_fit();
    INTER_POS.shrink_to_fit();
    #pragma omp parallel for
    for(uint32 i = 0; i < P_SIZE; ++i){
      uint i4 = i*4;
      uint i3 = i*3;
      POINT_DATA[i4] = points[i3];
      POINT_DATA[i4+1] = points[i3+1];
      POINT_DATA[i4+2] = points[i3+2];
      POINT_DATA[i4+3] = 0;
    }

    //cerr<<"Build log: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";


    cl::Buffer coord_data_cl(context,
                             CL_MEM_READ_ONLY,
                             sizeof(float)*COORD_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_2: " << err << "\n";
    cl::Buffer point_data_cl(context,
                             CL_MEM_READ_ONLY,
                             sizeof(float)*POINT_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_3: " << err << "\n";
    cl::Buffer int_data_cl(context,
                           CL_MEM_READ_ONLY,
                           sizeof(uint)*INT_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_4: " << err << "\n";
    cl::Buffer inter_pos_cl(context,
                            CL_MEM_READ_ONLY,
                            sizeof(uint)*INTER_POS.size(),NULL,&err);
//   cerr << "SimpleSolver_5: " << err << "\n";
    cl::Buffer area_count_cl(context,
                            CL_MEM_READ_WRITE,
                            sizeof(uint)*AREA_COUNT.size(),NULL,&err);
//   cerr << "area_count: " << err << "\n";
    cl::CommandQueue queue(context,
                           default_device);

    err = queue.enqueueWriteBuffer(coord_data_cl,
                                  CL_FALSE,
                                  0,
                                  sizeof(float)*COORD_DATA.size(),
                                  COORD_DATA.data());
//   cerr << "SimpleSolver_7: " << err << "\n";
    err = queue.enqueueWriteBuffer(point_data_cl,
                           CL_FALSE,
                           0,
                           sizeof(float)*POINT_DATA.size(),
                           POINT_DATA.data());
//   cerr << "SimpleSolver_8: " << err << "\n";
    err = queue.enqueueWriteBuffer(int_data_cl,
                           CL_FALSE,
                           0,
                           sizeof(uint)*INT_DATA.size(),
                           INT_DATA.data());
//   cerr << "SimpleSolver_9: " << err << "\n";
    err = queue.enqueueWriteBuffer(inter_pos_cl,
                           CL_FALSE,
                           0,
                           sizeof(uint)*INTER_POS.size(),
                           INTER_POS.data());
//   cerr << "SimpleSolver_10: " << err << "\n";

    err = queue.enqueueWriteBuffer(area_count_cl,
                                 CL_FALSE,
                                 0,
                                 sizeof(uint)*AREA_COUNT.size(),
                                 AREA_COUNT.data());
//   cerr << "SimpleSolver_12: " << err << "\n";
    cl::Kernel kernel = cl::Kernel(program,"simple_solver");
    kernel.setArg(0,coord_data_cl);
    kernel.setArg(1,point_data_cl);
    kernel.setArg(2,int_data_cl);
    kernel.setArg(3,inter_pos_cl);
    kernel.setArg(4,area_count_cl);
    queue.finish();
    queue.enqueueNDRangeKernel(kernel,
                              cl::NullRange,
                              cl::NDRange(pdb.size()*P_SIZE),
                              cl::NullRange);
    queue.enqueueReadBuffer(area_count_cl,
                            CL_TRUE,
                            0,
                            sizeof(uint)*AREA_COUNT.size(),
                            AREA_COUNT.data());
    
    #pragma omp parallel for
    for(uint i = 0; i < pdb.size(); ++i){
      //cout << pdb[i].ID << "\t" << pdb[i].STRUCT_TYPE <<"\t" << pdb[i].AREA << "\t" << AREA_COUNT[i] << "\t" << P_SIZE <<"\n";
      pdb[i].SASA = pdb[i].AREA*AREA_COUNT[i]/P_SIZE;
    }
  //end of FP64 branch
  }
  else {
    vector<float> COORD_DATA(pdb.size()*4);
    vector<float> POINT_DATA(P_SIZE*4);
    vector<uint> AREA_COUNT(pdb.size(),P_SIZE);
    vector<uint> INT_DATA;
    vector<uint> INTER_POS;
    //vector<char> AREAP(pdb.size()*P_SIZE,0);
    #pragma omp parallel for
    for (uint32 i = 0; i < pdb.size(); ++i){
      int i4 = i*4;
      COORD_DATA[i4] =     (float)pdb[i].RADIUS;
      COORD_DATA[i4 + 1] = (float)pdb[i].COORDS[0];
      COORD_DATA[i4 + 2] = (float)pdb[i].COORDS[1];
      COORD_DATA[i4 + 3] = (float)pdb[i].COORDS[2];
    }

    INTER_POS.push_back(0);
    for(auto& atom : pdb){
      INTER_POS.push_back(INTER_POS.back()+atom.INTERACTION_P.size());
      for(auto& p : atom.INTERACTION_P) {
        INT_DATA.push_back(p);
      }
    }
    INT_DATA.shrink_to_fit();
    INTER_POS.shrink_to_fit();
    #pragma omp parallel for
    for(uint32 i = 0; i < P_SIZE; ++i){
      uint i4 = i*4;
      uint i3 = i*3;
      POINT_DATA[i4] = (float)points[i3];
      POINT_DATA[i4+1] = (float)points[i3+1];
      POINT_DATA[i4+2] = (float)points[i3+2];
      POINT_DATA[i4+3] = (float)0;
    }
  
    //cerr<<"Build log: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";


    cl::Buffer coord_data_cl(context,
                             CL_MEM_READ_ONLY,
                             sizeof(float)*COORD_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_2: " << err << "\n";
    cl::Buffer point_data_cl(context,
                             CL_MEM_READ_ONLY,
                             sizeof(float)*POINT_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_3: " << err << "\n";
    cl::Buffer int_data_cl(context,
                           CL_MEM_READ_ONLY,
                           sizeof(uint)*INT_DATA.size(),NULL,&err);
  //cerr << "SimpleSolver_4: " << err << "\n";
    cl::Buffer inter_pos_cl(context,
                            CL_MEM_READ_ONLY,
                            sizeof(uint)*INTER_POS.size(),NULL,&err);
//   cerr << "SimpleSolver_5: " << err << "\n";
    cl::Buffer area_count_cl(context,
                             CL_MEM_READ_WRITE,
                             sizeof(uint)*AREA_COUNT.size(),NULL,&err);
//   cerr << "area_count: " << err << "\n";
    cl::CommandQueue queue(context,
                         default_device);

    err = queue.enqueueWriteBuffer(coord_data_cl,
                                   CL_FALSE,
                                   0,
                                   sizeof(float)*COORD_DATA.size(),
                                   COORD_DATA.data());
//   cerr << "SimpleSolver_7: " << err << "\n";
    err = queue.enqueueWriteBuffer(point_data_cl,
                                   CL_FALSE,
                                   0,
                                   sizeof(float)*POINT_DATA.size(),
                                   POINT_DATA.data());
//   cerr << "SimpleSolver_8: " << err << "\n";
    err = queue.enqueueWriteBuffer(int_data_cl,
                                   CL_FALSE,
                                   0,
                                   sizeof(uint)*INT_DATA.size(),
                                   INT_DATA.data());
//   cerr << "SimpleSolver_9: " << err << "\n";
    err = queue.enqueueWriteBuffer(inter_pos_cl,
                                   CL_FALSE,
                                   0,
                                   sizeof(uint)*INTER_POS.size(),
                                   INTER_POS.data());
//   cerr << "SimpleSolver_10: " << err << "\n";

    err = queue.enqueueWriteBuffer(area_count_cl,
                                   CL_FALSE,
                                   0,
                                   sizeof(uint)*AREA_COUNT.size(),
                                   AREA_COUNT.data());
//   cerr << "SimpleSolver_12: " << err << "\n";
    cl::Kernel kernel = cl::Kernel(program,"simple_solver");
    kernel.setArg(0,coord_data_cl);
    kernel.setArg(1,point_data_cl);
    kernel.setArg(2,int_data_cl);
    kernel.setArg(3,inter_pos_cl);
    kernel.setArg(4,area_count_cl);
    queue.finish();
    queue.enqueueNDRangeKernel(kernel,
                              cl::NullRange,
                              cl::NDRange(pdb.size()*P_SIZE),
                              cl::NullRange);
    queue.enqueueReadBuffer(area_count_cl,
                            CL_TRUE,
                            0,
                            sizeof(uint)*AREA_COUNT.size(),
                            AREA_COUNT.data());
    
    #pragma omp parallel for
    for(uint i = 0; i < pdb.size(); ++i){
      pdb[i].SASA = pdb[i].AREA*AREA_COUNT[i]/P_SIZE;
    }
  }

  //chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  //auto time_span = chrono::duration_cast<chrono::duration<float>>(t2 - t1);
  //cerr << "SimpleSolverCL took " << time_span.count() << " seconds.\n";*/
}
