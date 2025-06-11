 #pragma once 
 
 // Call helpers
   template <typename O>
   void GpuParallelizationHelper::gpuAtomCall(O op) {
      // Assert that O is derived from Atom
      (void)static_cast<Atom*>((O*)0);

      // Setup size
      op.NM = N * M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N * M);
      atom_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void GpuParallelizationHelper::gpuSiteCall(O op) {
      // Assert that O is derived from Site
      (void)static_cast<Site*>((O*)0);

      // Setup size
      op.N = N;
      op.NH = NH;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N);
      site_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void GpuParallelizationHelper::gpuAtomSiteCall(O op) {
      // Assert that O is derived from AtomSite
      (void)static_cast<AtomSite*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim2d(&block, &grid, N, M);
      atom_site_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void GpuParallelizationHelper::gpuAtomSiteEnsembleCall(O op) {
      // Assert that O is derived from AtomSiteEnsemble
      (void)static_cast<AtomSiteEnsemble*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;
      op.NH = NH;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim2d(&block, &grid, N, M);
      atom_site_ensemble_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void GpuParallelizationHelper::gpuElementAxisSiteEnsembleCall(O op) {
      // Assert that O is derived from ElementAxisSiteEnsemble
      (void)static_cast<ElementAxisSiteEnsemble*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;
      op.NH = NH;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim3d(&block, &grid, 3, N, M);
      element_axis_site_ensemble_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   // Call helpers
   template <typename O>
   void GpuParallelizationHelper::gpuElementCall(O op) {
      // Assert that O is derived from Element
      (void)static_cast<Element*>((O*)0);

      // Setup size
      op.NM3 = N * M * 3;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N * M * 3);
      element_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
      // gpuErrchk(cudaPeekAtLastError());
      // gpuErrchk(cudaDeviceSynchronize());
   }

