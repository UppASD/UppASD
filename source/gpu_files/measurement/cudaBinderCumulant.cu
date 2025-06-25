#include "cudaBinderCumulant.cuh"
#include "fortranData.hpp"
#include "cudaParallelizationHelper.hpp"

// should spawn one thread only
__global__ void naiveBinderCumulantNoEnergy_kernel(const CudaTensor<real, 3> emomM,
                                            real temp,
                                            real mub,
                                            real k_bolt,
                                            BinderCumulantData* d);


BinderCumulant::BinderCumulant(const CudaTensor<real, 3>& emomM,
                               dim3 cudaThreads)
: emomM(emomM)
, cumu_buff_size(*FortranData::cumu_buff)
, N(*FortranData::Natom)
, M(*FortranData::Mensemble)
, cumu_step(*FortranData::cumu_step)
, temp(*FortranData::temperature)
, plotenergy(static_cast<bool>(*FortranData::plotenergy))
, threads(cudaThreads)
, blocks( (N + threads.x - 1) / threads.x, (M + threads.y - 1) / threads.y )
{
    cumu_buff.Allocate(cumu_buff_size);
    cumu_buff.zeros();
}


BinderCumulant::~BinderCumulant()
{
    cumu_buff.Free();
}


void BinderCumulant::measure(std::size_t mstep)
{
    if ((mstep /*-1 or not?? it differs in fortran code from mavg*/) % cumu_step != 0)
        return;


    cudaStream_t workStream = CudaParallelizationHelper::def.getWorkStream();
    if (plotenergy)
    {
        assert(false); // not yet implemented
    }
    else
    {
        naiveBinderCumulantNoEnergy_kernel<<<1, 1, 0, workStream>>>(
            emomM,
            temp,
            *FortranData::mub,
            *FortranData::k_bolt,
            cumu_buff.data() + cumu_buff_count
        );
    }
    cudaDeviceSynchronize(); // for printing without delay
    ++cumu_buff_count;
    if (cumu_buff_count >= cumu_buff_size)
    {
        // TODO: save to Fortran
        cumu_buff_count = 0;
    }
}


void BinderCumulant::flushMeasurements(std::size_t mstep)
{

}


__global__ void naiveBinderCumulantNoEnergy_kernel(const CudaTensor<real, 3> emomM,
                                                  real temp,
                                                  real mub,
                                                  real k_bolt,
                                                  BinderCumulantData* d)
{
    const uint tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= 1)
        return;

    const uint atoms = emomM.extent(1);
    const uint ensembles = emomM.extent(2);

    for (uint k = 0; k < ensembles; ++k)
    {
        real m[3] = {0};

        for (uint i = 0; i < atoms; ++i)
        {
            m[0] += emomM(0, i, k);
            m[1] += emomM(1, i, k);
            m[2] += emomM(2, i, k);
        }

        const real avrgme = norm(3, m) / atoms;
        const real avrgm2 = pow(avrgme, 2);
        const real avrgm4 = pow(avrgm2, 2);
        d->cumuw += real(1);

        d->avrgmcum = (d->avrgmcum * d->cumutotw + avrgme * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm2cum = (d->avrgm2cum * d->cumutotw + avrgm2 * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm4cum = (d->avrgm4cum * d->cumutotw + avrgm4 * d->cumuw ) / (d->cumutotw + d->cumuw);

        if (d->avrgm2cum > 0.0)
            d->binderc = real(1) - (d->avrgm4cum / real(3) / pow(d->avrgm2cum, 2));

        // For T=0, not the proper susceptibility
        const real extra = (temp == 0.0)? 1.0 : k_bolt * temp;
        d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / (k_bolt * extra);

        d->Navrgcum += 1;
        d->cumutotw += d->cumuw;
    }



    printf("%d\t%e\t%e\t%e\t%e\t%e\t\n",
           d->Navrgcum / ensembles,
           d->avrgmcum,
           d->avrgm2cum,
           d->avrgm4cum,
           d->binderc,
           d->pmsusc
    );

}