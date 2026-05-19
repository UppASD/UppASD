#include "gpuMeasurement.hpp"
#include "c_helper.h"
#include "stopwatchPool.hpp"
#include "gpuParallelizationHelper.hpp"
#include <iostream>
#include "gpu_wrappers.h"

#include "measurementQueue.hpp"
using ParallelizationHelper = GpuParallelizationHelper;
namespace mm = kernels::measurement;

namespace {
inline bool valid_ptr(const void* p) { return p != nullptr; }
}

GpuMeasurement::GpuMeasurement(const deviceLattice& gpuLattice,
                                const deviceEnergies& gpuEnergies,
                                 Tensor<real, 3>& f_emomM, 
                                 Tensor<real, 3>& f_emom,
                                 Tensor<real, 2>& f_mmom,
                                 Tensor<real, 3>& f_beff,
                                 MeasurementQueue& mq,
                                 bool p_do_jtensor,
                                 bool alwaysCopy
                               )
: gpuLattice(gpuLattice)
, gpuEnergies(gpuEnergies)
, do_gpu_measurements(*FortranData::do_cuda_measurements)
, mqueue(mq)
, N(gpuLattice.emomM.extent(1))
, M(gpuLattice.emomM.extent(2))
, NX(0) // TODO: dont hardcode these values, needs to be imported from Fortran
, NY(0)
, NZ(0)
, NT(0)
, NA(*FortranData::NA)
, Nchmax(*FortranData::Nchmax)
, Natom_full(*FortranData::Natom_full)
//, NX(128) // TODO: dont hardcode these values, needs to be imported from Fortran
//, NY(128)
//, NZ(1)
//, NT(1)
, do_jtensor(p_do_jtensor)
, measurementWriter(do_jtensor)
, cpuMeas(gpuLattice.emomM, gpuLattice.emom, gpuLattice.mmom, gpuLattice.beff, f_emomM, f_emom, f_mmom, f_beff, mq)
, workStream( ParallelizationHelperInstance.getWorkStream() )
, stopwatch(GlobalStopwatchPool::get("Gpu measurement"))
, do_avrg(*FortranData::do_avrg == 'Y')
, do_avrg_proj(*FortranData::do_avrg_proj)
, do_avrg_projch(*FortranData::do_avrg_projch == 'Y')
, do_cumu_proj(*FortranData::do_cumu_proj)
, do_cumu_projch(0)
, do_ralloy(*FortranData::do_ralloy)
, mavg_kernel_threads(256)
, mavg_kernel_blocks(mm::ceil_div(M, mavg_kernel_threads.x))
, do_cumu(*FortranData::do_cumu == 'Y')
, cumu_kernel_threads(32)
, cumu_kernel_blocks(mm::ceil_div(M, cumu_kernel_threads.x))
, cumu_ene_maxBlocks(1024)
, cumu_ene_maxThreads(256)
, cumu_ene_kernel_threads(cumu_ene_maxThreads)
, cumu_ene_kernel_blocks(std::min(((M + cumu_ene_maxThreads - 1) / cumu_ene_maxThreads), cumu_ene_maxBlocks), 3, 1)
, sumOverAtoms_kernel_threads(256, 1)
, sumOverAtoms_kernel_blocks(mm::ceil_div(N, sumOverAtoms_kernel_threads.x), 3 * M)
, do_autocorr((*FortranData::do_autocorr=='Y') ? true : false)
, nspinwait(*FortranData::nspinwait)
, n0spinwait(0)
, ac_buff(*FortranData::ac_buff)
, ac_step(*FortranData::ac_step)
, sw_next(0)
, ac_maxThreads(256)
, ac_maxBlocks(1024)
, ac_count(0)
, do_skyno(skyrmionMethodFromFlag(*FortranData::do_skyno))
, nsimp(2 * NX * NY * NZ * NT)
, skyno_kernel_threads(128, 1)
, skyno_kernel_blocks(skyrmionKernelNumBlocks(do_skyno, N, M, nsimp, skyno_kernel_threads.x))
, do_ene(*FortranData::do_ene)
, ene_step(*FortranData::ene_step)
, ene_buff(*FortranData::ene_buff)
, ene_types(6)
, mub(9.274009994e-24)
, mry(2.179872325e-21)
, fcinv(mub / mry)
{
    
    isAllocated = false;
    asitealloc = false;

    if (do_avrg)
    {
        assert(*FortranData::avrg_step > 0 && *FortranData::avrg_buff > 0);

        mavg_buff_gpu.Allocate(*FortranData::avrg_buff);
        mavg_buff_gpu.zeros();

        mavg_buff_cpu.AllocateHost(*FortranData::avrg_buff);
        mavg_buff_cpu.zeros();

        mavg_partial_buff.Allocate(mavg_kernel_blocks.x);

        mavg_iter.AllocateHost(*FortranData::avrg_buff);
        mavg_iter.zeros();
        if((do_avrg_proj=='Y')){
            mavg_proj_kernel_blocks = {1,1,1};
            mavg_proj_buff_gpu.Allocate(NT, *FortranData::avrg_buff);
            mavg_proj_buff_gpu.zeros();
            mavg_proj_buff_cpu.AllocateHost(NT, *FortranData::avrg_buff);
            mavg_proj_buff_cpu.zeros();
            mavg_proj_partial_buff.Allocate(mavg_proj_kernel_blocks.x, NT);
            mavg_proj_partial_buff.zeros();
        }
        else if((do_avrg_proj=='A')){
            mavg_proj_kernel_blocks = {1,1,1};
            mavg_proj_buff_gpu.Allocate(NA, *FortranData::avrg_buff);
            mavg_proj_buff_gpu.zeros();
            mavg_proj_buff_cpu.AllocateHost(NA, *FortranData::avrg_buff);
            mavg_proj_buff_cpu.zeros();
            mavg_proj_partial_buff.Allocate(mavg_proj_kernel_blocks.x, NA);
            mavg_proj_partial_buff.zeros();

        }
        if((do_avrg_projch=='Y')){
            mavg_projch_buff_gpu.Allocate(Nchmax, *FortranData::avrg_buff);
            mavg_projch_buff_gpu.zeros();
            mavg_projch_buff_cpu.AllocateHost(Nchmax, *FortranData::avrg_buff);
            mavg_projch_buff_cpu.zeros();
            mavg_projch_partial_buff.Allocate(mavg_projch_kernel_blocks.x, Nchmax);
            mavg_projch_partial_buff.zeros();
        }

        std::cout << "AverageMagnetization observable added" << std::endl;
    }

    if (do_cumu)
    {
        assert(*FortranData::cumu_step > 0 && *FortranData::cumu_buff > 0);

        cumu_buff_gpu.Allocate(1);
        cumu_buff_gpu.zeros();

        cumu_buff_cpu.AllocateHost(1);
        cumu_buff_cpu.zeros();

        if (do_ene==0) cumu_partial_buff.Allocate(cumu_kernel_blocks.x);
        else cumu_ene_partial_buff.Allocate(cumu_ene_kernel_blocks.x);

        if(do_cumu_proj == 'Y'){
            cumu_proj_kernel_blocks = {1,1,1};

            cumu_proj_buff_gpu.Allocate(NT);
            cumu_proj_buff_gpu.zeros();

            cumu_proj_buff_cpu.AllocateHost(NT);
            cumu_proj_buff_cpu.zeros();

            cumu_proj_partial_buff.Allocate(cumu_proj_kernel_blocks.x, NT);
            cumu_proj_partial_buff.zeros();
        }
        else if(do_cumu_proj == 'A'){
            cumu_proj_kernel_blocks = {1,1,1};

            cumu_proj_buff_gpu.Allocate(NA);
            cumu_proj_buff_gpu.zeros();

            cumu_proj_buff_cpu.AllocateHost(NA);
            cumu_proj_buff_cpu.zeros();

            cumu_proj_partial_buff.Allocate(cumu_proj_kernel_blocks.x, NA);
            cumu_proj_partial_buff.zeros();

        }
        
        if(do_cumu_projch){

            cumu_projch_kernel_blocks = {1,1,1};

            cumu_projch_buff_gpu.Allocate(NT);
            cumu_projch_buff_gpu.zeros();

            cumu_projch_buff_cpu.AllocateHost(NT);
            cumu_projch_buff_cpu.zeros();

            cumu_projch_partial_buff.Allocate(cumu_projch_kernel_blocks.x, NT);
            cumu_projch_partial_buff.zeros();
        }

        std::cout << "BinderCumulant observable added" << std::endl;
    }

    if (do_avrg || do_cumu)
    {
        emomMEnsembleSums.Allocate(3, M);
        emomMEnsembleSums.zeros();

        emomMEnsembleSums_partial.Allocate(sumOverAtoms_kernel_blocks.x, 3, M);
        emomMEnsembleSums_partial.zeros();

        if((do_avrg_proj=='Y')||(do_cumu_proj == 'Y')){
            emomMEnsembleNTSums.Allocate(3, NT, M);
            emomMEnsembleNTSums.zeros();
            emomMEnsembleNTSums_partial.Allocate(sumOverAtoms_NT_kernel_blocks.x, 3, NT, M);
            emomMEnsembleNTSums_partial.zeros();
            atype_gpu.Allocate(N);
            atype_cpu.set(FortranData::atype, static_cast<long int>(N));
            atype_gpu.copy_sync(atype_cpu);

        }
        if((do_avrg_proj=='A')||(do_cumu_proj=='A')){
            emomMEnsembleNASums.Allocate(3, NA, M);
            emomMEnsembleNASums.zeros();
            emomMEnsembleNASums_partial.Allocate(sumOverAtoms_NA_kernel_blocks.x, 3, NA, M);
            emomMEnsembleNASums_partial.zeros();
            if((!asitealloc)&&(!do_ralloy) && valid_ptr(FortranData::asite_ch)){
                asite_ch_gpu.Allocate(N);
                asite_ch_cpu.set(FortranData::asite_ch, static_cast<long int>(N));
                asite_ch_gpu.copy_sync(asite_ch_cpu);
                asitealloc = true;
            }

        }
        if((do_avrg_projch=='Y')||(do_cumu_projch=='Y')){
            emomMEnsembleNCSums.Allocate(3, Nchmax, M);
            emomMEnsembleNCSums.zeros();
            emomMEnsembleNCSums_partial.Allocate(sumOverAtoms_NC_kernel_blocks.x, 3, Nchmax, M);
            emomMEnsembleNCSums_partial.zeros();
            achem_ch_gpu.Allocate(N);
            achem_ch_cpu.set(FortranData::achem_ch, static_cast<long int>(N));
            achem_ch_gpu.copy_sync(achem_ch_cpu);
            
        }
    }

    if (do_skyno == SkyrmionMethod::BruteForce || do_skyno == SkyrmionMethod::Triangulation)
    {
        assert(*FortranData::skyno_step > 0 && *FortranData::skyno_buff > 0);

        skyno_buff_gpu.Allocate(*FortranData::skyno_buff);
        skyno_buff_gpu.zeros();

        skyno_buff_cpu.AllocateHost(*FortranData::skyno_buff);
        skyno_buff_cpu.zeros();

        skyno_iter.AllocateHost(*FortranData::skyno_buff);
        skyno_iter.zeros();

        std::cout << "SkyrmionNumber observable added" << std::endl;
    }

    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        Tensor<real, 3> dxyz_vec_fortran(FortranData::dxyz_vec, 3, 26, N);
        dxyz_vec.Allocate(3, 26, N);
        dxyz_vec.copy_sync(dxyz_vec_fortran);

        Tensor<int, 2> dxyz_atom_fortran(FortranData::dxyz_atom, 26, N);
        dxyz_atom.Allocate(26, N);
        dxyz_atom.copy_sync(dxyz_atom_fortran);

        Tensor<int, 1> dxyz_list_fortran(FortranData::dxyz_list, N);
        dxyz_list.Allocate(N);
        dxyz_list.copy_sync(dxyz_list_fortran);

        grad_mom.Allocate(3, 3, N, M);
        grad_mom.zeros();

        skyno_partial_buff.Allocate(skyno_kernel_threads.x * skyno_kernel_threads.y);
    }

    if (do_skyno == SkyrmionMethod::Triangulation)
    {
        if (nsimp == 0)
            throw std::invalid_argument("NX, NY, NZ and NT needs to be set from Fortran. Currently only works with hard coded values.");

        simp.Allocate(3, nsimp);
        skyno_partial_buff.Allocate(skyno_kernel_blocks.x);

        const uint pairs = nsimp / 2;
        const uint threads = 256;
        const uint blocks = mm::ceil_div(pairs, threads);

        mm::delaunay_tri_tri<<<blocks, threads, 0, workStream>>>(NX, NY, NZ, NT, simp);
    }

    if (do_ene>0)
    {

        assert(*FortranData::ene_step > 0 && *FortranData::ene_buff > 0);

        ene_maxBlocks = 1024;
        ene_maxThreads = 256;
        ene_kernel_threads = ene_maxThreads;
        ene_kernel_blocks =  dim3(std::min(((M + ene_maxThreads - 1) / ene_maxThreads), ene_maxBlocks), ene_types, 1);
        energy_buff_gpu.Allocate(*FortranData::ene_buff);
        energy_buff_cpu.AllocateHost(*FortranData::ene_buff);
        energy_buff_gpu.zeros();
        energy_buff_cpu.zeros();
        energy_iter.AllocateHost(*FortranData::ene_buff);
        energy_iter.zeros();
        energy_partial_buff.Allocate(ene_kernel_blocks.x, ene_types);

        printf("WARNING: DO NOT USE BIG_GRID WITH GPU CALCULATIONS OF ENERGY\n");
    }
    if (do_autocorr){
        sw_threads = 256;
        sw_tasks = 3*N*M;
        sw_blocks = (sw_tasks + sw_threads - 1)/sw_threads;
        indxb_ac.set(FortranData::indxb_ac, ac_buff);
        //spinwait_cpu.set (FortranData::spinwait, 3, N, nspinwait);
        spinwaittable_cpu.set(FortranData::spinwaittable, nspinwait);
        spinwait_gpu.Allocate(3, N, M, nspinwait);
        spinwait_gpu.zeros();
        fill_spinwait<<<sw_blocks, sw_threads>>>(spinwait_gpu, gpuLattice.emom, sw_tasks, 0);
        
        sw_curIdx = 0;
        printf("\n do_ac = %c, swtt size = %i\n ", do_autocorr, spinwaittable_cpu.extent(0));
        sw_next = spinwaittable_cpu(0);
        sw_curr = 0;

        ac_tasksX = 3 * N * M;
        ac_threadsX = ac_maxThreads;
        ac_blocksX =  std::min(((ac_tasksX + ac_threadsX - 1) / ac_threadsX), ac_maxBlocks);
        ac_threads = {ac_threadsX, 1, 1};

        ac_block_gpu.Allocate(ac_blocksX, nspinwait);
        ac_block_gpu.zeros();

        autocorr_buff_gpu.Allocate(nspinwait, ac_buff);
        autocorr_buff_gpu.zeros();

        autocorr_buff_cpu.AllocateHost(nspinwait, ac_buff);
        autocorr_buff_cpu.zeros(); 
        
        indxb_ac.AllocateHost(ac_buff);
        indxb_ac.zeros();
    }


    isAllocated = true;
    stopwatch.add("constructor");
}

void GpuMeasurement::release(){
    if(!isAllocated)
        return;

    if (do_avrg)
    {
        mavg_buff_gpu.Free();
        mavg_buff_cpu.FreeHost();
        mavg_partial_buff.Free();
        mavg_iter.FreeHost();
        if((do_avrg_proj=='Y')||(do_avrg_proj=='A')){
            mavg_proj_buff_gpu.Free();
            mavg_proj_buff_cpu.FreeHost();
            mavg_proj_partial_buff.Free();
        }
        if((do_avrg_projch=='Y')){
            mavg_projch_buff_gpu.Free();
            mavg_projch_buff_cpu.FreeHost();
            mavg_projch_partial_buff.Free();
        }
    }

    if (do_cumu)
    {
        cumu_buff_gpu.Free();
        cumu_buff_cpu.FreeHost();
        if (do_ene==0) cumu_partial_buff.Free();
        else cumu_ene_partial_buff.Free();

        if((do_cumu_proj == 'Y')||(do_cumu_proj == 'A')){
            cumu_proj_buff_gpu.Free();
            cumu_proj_buff_cpu.FreeHost();
            cumu_proj_partial_buff.Free();
        }

        if(do_cumu_projch){
            cumu_projch_buff_gpu.Free();
            cumu_projch_buff_cpu.FreeHost();
            cumu_projch_partial_buff.Free();
        }

    }

    if (do_avrg || do_cumu)
    {
        emomMEnsembleSums.Free();
        emomMEnsembleSums_partial.Free();

        if((do_avrg_proj=='Y')||(do_cumu_proj == 'Y')){
            emomMEnsembleNTSums.Free();
            emomMEnsembleNTSums_partial.Free();
            atype_gpu.Free();
            if((asitealloc)&&(!do_ralloy)){
                asite_ch_gpu.Free();
                asitealloc = false;
            }
        }
        if((do_avrg_proj=='A')||(do_cumu_proj=='A')){
            emomMEnsembleNASums.Free();
            emomMEnsembleNASums_partial.Free();
            if((asitealloc)&&(!do_ralloy)){
                asite_ch_gpu.Free();
                asitealloc = false;
            }
        }
        if((do_avrg_projch=='Y')||(do_cumu_projch=='Y')){
            emomMEnsembleNCSums.Free();
            emomMEnsembleNCSums_partial.Free();
            achem_ch_gpu.Free();          
        }     
    }

    if (do_skyno == SkyrmionMethod::BruteForce || do_skyno == SkyrmionMethod::Triangulation)
    {
        skyno_buff_gpu.Free();
        skyno_buff_cpu.FreeHost();
        skyno_iter.FreeHost();
        skyno_partial_buff.Free();
    }

    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        dxyz_vec.Free();
        dxyz_atom.Free();
        dxyz_list.Free();
        grad_mom.Free();
    }

    if (do_skyno == SkyrmionMethod::Triangulation)
    {
        simp.Free();
    }

    if (do_ene>0)
    {
        energy_buff_gpu.Free();
        energy_partial_buff.Free();
        energy_buff_cpu.FreeHost();
        energy_iter.FreeHost();
    }
      
    if (do_autocorr)
    {
        spinwait_gpu.Free();
        ac_block_gpu.Free();
        autocorr_buff_gpu.Free();
        autocorr_buff_cpu.FreeHost();
        indxb_ac.FreeHost();
    }

    isAllocated = false;
}

GpuMeasurement::~GpuMeasurement()
{
    release();
}


void GpuMeasurement::measure(std::size_t mstep)
{
    cpuMeas.measure(mstep);

    --mstep; // this is because the simulation loop begins at 1 because of Fortran indexing

    const bool avrg = timeToMeasure(MeasurementType::AverageMagnetization, mstep);
    const bool cumu = timeToMeasure(MeasurementType::BinderCumulant, mstep);
    const bool energy = timeToMeasure(MeasurementType::Energy, mstep);
    const bool autocorr = timeToMeasure(MeasurementType::Autocorrelation, mstep);

    if (avrg || cumu)
    {
        calculateEmomMSum();
        stopwatch.add("sum reduction of emomM for shared use");
    }

    if (avrg)
    {
        measureAverageMagnetization(mstep);
        stopwatch.add("average magnetization");
    }

    if (energy)
    {
        measureEnergy(mstep);
        stopwatch.add("energy");
    }

    if (cumu)
    {
        measureBinderCumulant(mstep);
        stopwatch.add("binder cumulant");
    }

    if (timeToMeasure(MeasurementType::SkyrmionNumber, mstep))
    {
        measureSkyrmionNumber(mstep);
        stopwatch.add("skyrmion number");
    }
    if(autocorr){
        measureAutocorrelation(mstep);
        stopwatch.add("autocorelation");

    }//0 - exch, 1 - ani, 2 - dm, 3 - tensor, 4 - external, 5 - total

    

    if(GPU_DEVICE_SYNCHRONIZE() != GPU_SUCCESS) {
      release();
    }

}


void GpuMeasurement::flushMeasurements(std::size_t mstep)
{
    cpuMeas.measure(mstep);  

    --mstep; // this is because the simulation loop begins at 1 because of Fortran indexing

    if (do_avrg || do_cumu)
    {
        calculateEmomMSum();
        stopwatch.add("sum reduction of emomM for shared use");
    }

    if (do_avrg)
    {
        measureAverageMagnetization(mstep);
        stopwatch.add("average magnetization");

        if (mavg_count > 0)
            saveToFile(MeasurementType::AverageMagnetization);
    }

    if (do_gpu_measurements && (do_ene==1))
    {
        measureEnergy(mstep);
        stopwatch.add("energy");

        if (energy_count > 0)
            saveToFile(MeasurementType::Energy);
    }

    if (do_cumu)
    {
        measureBinderCumulant(mstep);
        stopwatch.add("binder cumulant");

        if (cumu_count > 0)
            saveToFile(MeasurementType::BinderCumulant);
    }

    if (do_skyno != SkyrmionMethod::None)
    {
        measureSkyrmionNumber(mstep);
        stopwatch.add("skyrmion number");

        if (skyno_count > 0)
            saveToFile(MeasurementType::SkyrmionNumber);
    }

    if (do_autocorr)
    {
        //measureAutocorrelation(mstep);
        //ac_count++;
        //stopwatch.add("autocorrelation");
        //saveToFile(MeasurementType::Autocorrelation);
    }

    cpuMeas.flushMeasurements(mstep + 1); 

    if(GPU_DEVICE_SYNCHRONIZE() != GPU_SUCCESS)
    {
        release();
    }
}


void GpuMeasurement::measureAverageMagnetization(std::size_t mstep)
{
    const size_t smem = mm::nwarps(mavg_kernel_threads) * sizeof(real);

    mm::averageMagnetization_partial<<<mavg_kernel_blocks, mavg_kernel_threads, smem, workStream>>>(
            emomMEnsembleSums, N, M, mavg_partial_buff.data()
    );


    mm::averageMagnetization_finalize<<<1, mavg_kernel_threads, smem, workStream>>>(
            mavg_partial_buff.data(), mavg_kernel_blocks.x, M, mavg_buff_gpu.data()[mavg_count]
    );

    fill_index(mavg_iter, mstep, mavg_count);
    mavg_count++;

    if (mavg_count >= *FortranData::avrg_buff)
    {
        saveToFile(MeasurementType::AverageMagnetization);
    }
}


void GpuMeasurement::measureBinderCumulant(std::size_t mstep)
{
    if ((do_ene==0))
    {
        const size_t smem = mm::nwarps(cumu_kernel_threads) * sizeof(real);

        mm::binderCumulantNoEnergy_partial<<<cumu_kernel_blocks, cumu_kernel_threads, smem, workStream>>>(
                emomMEnsembleSums, N, M, cumu_partial_buff.data()
        );

        mm::binderCumulantNoEnergy_finalize<<<1, cumu_kernel_threads, smem, workStream>>>(
                cumu_partial_buff.data(),
                cumu_kernel_blocks.x,
                N,
                M,
                *FortranData::temperature,
                *FortranData::mub,
                *FortranData::k_bolt,
                *cumu_buff_gpu.data()
        );
    }
    else
    {
        //throw std::invalid_argument("Not yet implemented.");

        mm::binderCumulantEnergy_partial<<<cumu_ene_kernel_blocks, cumu_ene_kernel_threads>>>(
                emomMEnsembleSums, gpuEnergies.energyM, N, M, cumu_ene_partial_buff.data()
        );

        mm::binderCumulantEnergy_finalize<<<1, cumu_ene_maxBlocks>>>(
                                                                        cumu_ene_partial_buff.data(),
                                                                        cumu_ene_kernel_blocks.x,
                                                                        N,
                                                                        M,
                                                                        *FortranData::temperature,
                                                                        *FortranData::mub,
                                                                        *FortranData::k_bolt,
                                                                        *FortranData::mry,
                                                                        *cumu_buff_gpu.data()
                                                                    );
    }


    if ((cumu_count++ % *FortranData::cumu_buff) == 0)
    {
        saveToFile(MeasurementType::BinderCumulant);
    }

}

void GpuMeasurement::measureSkyrmionNumber(std::size_t mstep)
{
    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        mm::grad_moments<<<skyno_kernel_blocks, skyno_kernel_threads, 0, workStream>>>(
                gpuLattice.emomM, dxyz_vec, dxyz_atom, dxyz_list, grad_mom
        );


        size_t smem = mm::nwarps(skyno_kernel_threads) * sizeof(real);
        mm::pontryagin_no_partial<<<skyno_kernel_blocks, skyno_kernel_threads, smem, workStream>>>(
                gpuLattice.emomM, grad_mom, skyno_partial_buff.data()
        );

        smem = skyno_kernel_threads.x * sizeof(real);
        mm::pontryagin_no_finalize<<<1, skyno_kernel_threads, smem, workStream>>>(
                skyno_partial_buff.data(), skyno_kernel_blocks.x, M, skyno_count + 1, skyno_buff_gpu.data()[skyno_count]
        );
    }

    else if (do_skyno == SkyrmionMethod::Triangulation)
    {
        size_t smem = mm::nwarps(skyno_kernel_threads) * sizeof(real);
        mm::pontryagin_tri_partial<<<skyno_kernel_blocks, skyno_kernel_threads, smem, workStream>>>(
                gpuLattice.emom, simp, skyno_partial_buff.data()
        );

        smem = skyno_kernel_threads.x * sizeof(real);
        mm::pontryagin_tri_finalize<<<1, skyno_kernel_threads, smem, workStream>>>(
                skyno_partial_buff.data(), skyno_kernel_blocks.x, M, skyno_count + 1, skyno_buff_gpu.data()[skyno_count]
        );
    }

    fill_index(skyno_iter, mstep, skyno_count);
    skyno_count++;

    if (skyno_count >= *FortranData::skyno_buff)
    {
        saveToFile(MeasurementType::SkyrmionNumber);
    }
}


void GpuMeasurement::measureEnergy(size_t mstep)
{
    // copy from gpuLattice to energy gpu buffer,
    // all the calculations for energy is done together with hamiltonian calculations
    //GPU_MEMCPY(energy_buff_gpu.data()+energy_count, gpuLattice.energy.data(),
           // 1 * sizeof(EnergyData), GPU_MEMCPY_DEVICE_TO_DEVICE);
           
           
    mm::averageEnergy_partial<<<ene_kernel_blocks, ene_kernel_threads>>>(gpuEnergies.energyM, M, energy_partial_buff.data());

    mm::averageEnergy_final<<<1, ene_maxBlocks>>>(
            energy_partial_buff.data(), ene_kernel_blocks.x, M, fcinv, energy_buff_gpu.data()[energy_count]);

   //printf("ene_step = %i, mstep = %i, ene_buff = %i, ene_ext = %i, ene_count = %i\n", 
    //ene_step, mstep, ene_buff,  energy_buff_gpu.extent(0), energy_count);
    
    fill_index(energy_iter, mstep, energy_count);
    energy_count++;

    if ((energy_count % *FortranData::ene_buff) == 0)
    {
        saveToFile(MeasurementType::Energy);
    }
}

void GpuMeasurement::measureAutocorrelation(std::size_t mstep)
{
 
   // printf("sw_next = %i, mstep = %i\n", sw_next, mstep);

    if (((mstep % ac_step) == 0)||(mstep == sw_curr)){
        ac_blocks = {ac_blocksX, sw_curIdx + 1, 1};
        real norm = 1/(static_cast<real>(gpuLattice.emom.extent(1)*gpuLattice.emom.extent(2)));
        //printf("cur swIdx = %i,  max sw = %i\n\n", sw_curIdx, nspinwait);
        //if (mstep == sw_curr) printf("HERE!\n");


        calc_autocorr_block<<<ac_blocks, ac_threads>>>(ac_block_gpu, spinwait_gpu, gpuLattice.emom);
        calc_autocorr_final<<<(sw_curIdx + 1), ac_maxBlocks>>>(ac_block_gpu, autocorr_buff_gpu, norm, ac_count, ac_blocksX);
    

        fill_index(indxb_ac, mstep + 1, ac_count);
        ac_count++;

        if (ac_count >= *FortranData::ac_buff)
        {
            saveToFile(MeasurementType::Autocorrelation);

        }
    }
   
}
   
void GpuMeasurement::updateAC(std::size_t mstep)
{
 
    if ((mstep == sw_next)&&(sw_curIdx < (spinwaittable_cpu.extent(0) - 1) )){
        printf("sw_next = %i, mstep = %i, swId = %i, \n\n", sw_next,  mstep, sw_curIdx);

        sw_curIdx++;
        fill_spinwait<<<sw_blocks, sw_threads>>>(spinwait_gpu, gpuLattice.emom, sw_tasks, sw_curIdx);
        sw_curr = sw_next;
        sw_next = spinwaittable_cpu(sw_curIdx); 
    }
   
}
   


void GpuMeasurement::calculateEmomMSum()
{
    size_t smem = mm::nwarps(sumOverAtoms_kernel_threads) * sizeof(real);
    mm::sumOverAtoms_partial<<<sumOverAtoms_kernel_blocks, sumOverAtoms_kernel_threads, smem, workStream>>>(
            gpuLattice.emomM,
            emomMEnsembleSums_partial
    );

    smem = sumOverAtoms_kernel_blocks.x * sizeof(real);
    const dim3 threads = 256;
    const dim3 grid = 3 * M;
    mm::sumOverAtoms_finalize<<<grid, threads, smem, workStream>>>(
            emomMEnsembleSums_partial,
            sumOverAtoms_kernel_blocks.x,
            emomMEnsembleSums
    );
}


void GpuMeasurement::saveToFile(MeasurementType mtype)
{
    switch (mtype)
    {
    case MeasurementType::AverageMagnetization:
        mavg_buff_cpu.copy_sync(mavg_buff_gpu);

        measurementWriter.write(
                mavg_iter.data(),
                mavg_buff_cpu.data(),
                mavg_count
        );

        mavg_buff_gpu.zeros();
        mavg_buff_cpu.zeros();
        mavg_count = 0;
        mavg_iter.zeros();
    break;

    case MeasurementType::BinderCumulant:
        cumu_buff_cpu.copy_sync(cumu_buff_gpu);

        measurementWriter.write(
                &cumu_count, // TODO: in fortran the equiv to cumu_count is printed, but should it not be mstep?
                cumu_buff_cpu.data(),
                1
        );
    break;

    case MeasurementType::SkyrmionNumber:
        skyno_buff_cpu.copy_sync(skyno_buff_gpu);

        measurementWriter.write(
                skyno_iter.data(),
                skyno_buff_cpu.data(),
                skyno_count
        );

        skyno_buff_gpu.zeros();
        skyno_buff_cpu.zeros();
        skyno_iter.zeros();
        skyno_count = 0;
    break;

    case MeasurementType::Energy:
        energy_buff_cpu.copy_sync(energy_buff_gpu);

        measurementWriter.write(
            energy_iter.data(),
            energy_buff_cpu.data(),
            energy_count
        );
        measurementWriter.writeEnergyStd(
            energy_iter.data(),
            energy_buff_cpu.data(),
            energy_count
        );

        energy_buff_gpu.zeros();
        energy_buff_cpu.zeros();
        energy_iter.zeros();
        energy_count = 0;
        //printf("HERE!, ene count = %i", energy_count);
    
    break;
        
    case MeasurementType::Autocorrelation:
        autocorr_buff_cpu.copy_sync(autocorr_buff_gpu);

        for (size_t i = 0; i < ac_buff; ++i)
        {
            AutocorrelationData row{
                .values = autocorr_buff_cpu.data() + i * nspinwait,
                .size   = nspinwait
            };
            
            measurementWriter.write(&indxb_ac[i], &row, 1);

        }
        autocorr_buff_gpu.zeros();
        autocorr_buff_cpu.zeros();
        indxb_ac.zeros();
        ac_count = 0;
  break;

  default:
      throw std::invalid_argument("Not yet implemented.");

    }
}


bool GpuMeasurement::timeToMeasure(MeasurementType mtype, size_t mstep) const
{
    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            return do_avrg && ((mstep % *FortranData::avrg_step) == 0);

        case MeasurementType::BinderCumulant:
            return do_cumu && ((mstep % *FortranData::cumu_step) == 0);

        case MeasurementType::SkyrmionNumber:
            return do_skyno != SkyrmionMethod::None && ((mstep % *FortranData::skyno_step) == 0);
        
        case MeasurementType::Autocorrelation: {
            return (do_autocorr=='Y') && (((mstep % ac_step) == 0)||(mstep== sw_curr));//TODO
        }

        case MeasurementType::Energy:
            return do_gpu_measurements && (do_ene>0) && ((mstep % *FortranData::ene_step) == 0);

        default:
            throw std::invalid_argument("Not yet implemented.");
    }
}


dim3 GpuMeasurement::skyrmionKernelNumBlocks(SkyrmionMethod method, uint N, uint M, uint nsimp, uint kernel_threads)
{
    switch (method)
    {
        case SkyrmionMethod::BruteForce: return (mm::ceil_div(N, kernel_threads), M);
        case SkyrmionMethod::Triangulation: return mm::ceil_div(nsimp, kernel_threads);
        default: return 0;
    }
}

template<typename T>
inline void GpuMeasurement::fill_index(Vector<real>& iter, T step, size_t& count)
{
    if (*FortranData::real_time_measure == 'Y') {
        iter(count) = static_cast<real>(step) * (*FortranData::delta_t);
    }
    else {
        iter(count) = static_cast<real>(step);
    }
}
