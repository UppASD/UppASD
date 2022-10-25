""" @package ASDPlots2D
Set of classes to perform the matplotlib plots of the UppASD data.

It contains two major classes:
    - Line/scatter plots for 2D and 3D data.
    - Spin-Spin correlation plots which require using colormaps.

Author
----------
Jonathan Chico
"""
################################################################################
# @brief Class defining the needed structures for line and scatter plots.
# @details This routine is used to plot 2D and 3D plots in matplotlib from
# \c UppASD data.
#
# @author Jonathan Chico
################################################################################
class Abstract2DPlot():
    def __init__(self):
        Abstract2DPlot.font_size=12
        Abstract2DPlot.markersize=0
        Abstract2DPlot.linewidth=3
        Abstract2DPlot.xgrid=False
        Abstract2DPlot.ygrid=False
        Abstract2DPlot.amsgrid=False
        return
    ############################################################################
    # @brief Function to plot line and scatter plots.
    # @details This function takes an array of arrays and plot them sequentially
    # the scatter points are located on top of solid lines, with colors been given
    # by a Paired colormap whose size is determined by the size of the array.
    # @author Jonathan Chico
    ############################################################################
    def LinePlot(self,axis,data_x,data_y,labels,ax_label,**kwargs):
        from matplotlib import cm as cm
        import numpy as np

        tick_labels = kwargs.get('tick_labels',None)
        tick_idx = kwargs.get('tick_idx',None)

        # AB -> add figure properties
        axis.cla()
        colors=cm.Paired(np.linspace(0,1,len(data_x)+2))
        for ii in range(0,len(data_x)):
            axis.plot(data_x[ii],data_y[ii],lw=self.linewidth,c=colors[ii],label=labels[ii],zorder=-1,markersize=self.markersize,marker='o')
            #axis.scatter(data_x[ii],data_y[ii],color=colors[ii],alpha=0.75, s=150,lw=1.00, edgecolor='black')

        axis.xaxis.grid(visible=self.xgrid,which='major')
        axis.yaxis.grid(visible=self.ygrid,which='major')
        axis.set_xlabel(ax_label[0],fontsize=Abstract2DPlot.font_size)
        axis.set_ylabel(ax_label[1],fontsize=Abstract2DPlot.font_size)
        axis.tick_params(axis='x', colors='black',which='minor',labelsize=Abstract2DPlot.font_size,width=2)
        axis.tick_params(axis='y', colors='black',which='minor',labelsize=Abstract2DPlot.font_size,width=2)
        axis.legend(fontsize=Abstract2DPlot.font_size)
        for ax in ['top','bottom','left','right']:
            axis.spines[ax].set_linewidth(3)

        if  self.amsgrid and tick_idx and tick_labels:
            axis.set_xticks(tick_idx,tick_labels)
            axis.xaxis.grid(visible=self.amsgrid,which='major',color='k')
            axis.autoscale(enable=True, axis='x', tight=True)
            axis.set_ylim([0,None])
        else:
            axis.xaxis.grid(visible=self.amsgrid)
            axis.autoscale(enable=True, axis='x', tight=False)

        return
    ############################################################################
    # @brief Function to plot a 3D representation of a magnetic moment trajectory.
    # @details This function takes an array of arrays and plot them sequentially
    # with colors been given by a Paired colormap whose size is determined by the
    # size of the array. The moments are constrained to the unit sphere.
    # @author Jonathan Chico
    ############################################################################
    def TrajPlot(self,axis,traj_data_x,traj_data_y,traj_data_z,traj_label):
        from matplotlib import cm as cm
        import numpy as np
        axis.cla()
        axis.set_xlim([-1,1])
        axis.set_ylim([-1,1])
        axis.set_zlim([-1,1])
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        axis.plot_wireframe(x, y, z, color="gray", linewidth=0.1)
        axis.set_xticks([-1,0,1])
        axis.set_yticks([-1,0,1])
        axis.set_zticks([-1,0,1])
        colors=cm.Paired(np.linspace(0,1,len(traj_data_x)+2))
        for ii in range(len(traj_data_x)):
            axis.plot(traj_data_x[ii],traj_data_y[ii],traj_data_z[ii], c=colors[ii],lw=3,label=traj_label[ii])
        axis.legend(fontsize=Abstract2DPlot.font_size)
        axis.tick_params(axis='x', colors='black',labelsize=Abstract2DPlot.font_size,width=2)
        axis.tick_params(axis='y', colors='black',labelsize=Abstract2DPlot.font_size,width=2)
        axis.tick_params(axis='z', colors='black',labelsize=Abstract2DPlot.font_size,width=2)
        axis.set_xlabel(r'M$_x$',fontsize=Abstract2DPlot.font_size)
        axis.set_ylabel(r'M$_y$',fontsize=Abstract2DPlot.font_size)
        axis.set_zlabel(r'M$_z$',fontsize=Abstract2DPlot.font_size)
        return
################################################################################
# @brief Class defining the needed structures for the spin-spin correlation plots.
#
# @ details Due to the fact that these make heavy use of colormap structures they are kept
# separated from the rest of the plotting routines.
#
# @author Jonathan Chico
################################################################################
class Correlation_Plots():
    ############################################################################
    # @brief Constructor of the class
    # @details It contains the different colormaps used for the \f$\mathbf{S}(\mathbf{q},\omega)\f$ plotting
    # @author Jonathan Chico
    ############################################################################
    def __init__(self):
        from matplotlib import cm as cm
        Correlation_Plots.font_size=12
        Correlation_Plots.cmap=[cm.coolwarm,cm.Spectral,cm.inferno]
        Correlation_Plots.sigma_w=0.0
        Correlation_Plots.w_max=1.0
        Correlation_Plots.w_min=1.0
        Correlation_Plots.grid = False
        return
    ############################################################################
    # @brief Create a Gaussian kernel (normalized) for future convolutions
    # @author Anders Bergman
    ############################################################################
    def gaussian(self,M,std):
        import numpy as np
        x=np.linspace(-self.w_max/2.0,self.w_max/2.0,M)
        y=(1.0)/(std*np.sqrt(2.0*np.pi))*np.exp(-0.5*x*x/std**2)
        # For consistency with scipy.signal - skip normalization as per below
        # y_0=np.exp(-0.5*x*x/std**2)
        return y

    ############################################################################
    # @brief Convolute a function with a kernel
    # @author Anders Bergman
    ############################################################################
    def convolve(self,func,kern):
        import numpy as np
        n=func.shape[0]+kern.shape[0]
        F=np.fft.fft(func,n=n)
        G=np.fft.fft(kern,n=n)
        FG=F*G
        zpos=int((kern.shape[0]-1)/2.0)
        return np.real(np.fft.ifft(FG))[zpos:zpos+func.shape[0]]

    def sqw_convolve(self,sqwa):
        import numpy as np
        #-----------------------------------------------------------------------
        # Perform a convolution with a windowing function for each q-point
        #-----------------------------------------------------------------------
        ed=sqwa.shape[0]
        qd=sqwa.shape[1]
        gauss=self.gaussian(ed,self.sigma_w)
        sqw_in=np.copy(sqwa)
        for iq in range(0,qd):
            sqw_in[:,iq]=self.convolve(sqw_in[:,iq],gauss)
        #-----------------------------------------------------------------------
        # Find the peaks and normalize the data
        #-----------------------------------------------------------------------
        sqw_peaks=np.argmax(sqw_in,axis=0)
        normMat=np.diag(1.0/np.amax(sqw_in,axis=0))
        sqw_out=np.matmul(sqw_in,normMat)
        return sqw_out

    ############################################################################
    # @brief Function to plot the spin-spin correlation function obtained from the
    # sqw.simid.out file.
    #
    # @author Jonathan Chico
    ############################################################################
    def Sqw_Plot(self,axis,sqw_data,proj,sqw_labels,col_indx,ax_limits,q_labels,q_idx):
        axis.cla()
        import numpy as np
        # Energy maximim 
        self.w_max=ax_limits[3]
        self.w_min=self.w_max/sqw_data[proj].shape[0]
        # Setting the Gaussian smearing to minimum relevant value (0 is not good..)
        self.sigma_w=np.maximum(self.sigma_w,self.w_max/sqw_data[proj].shape[0])
        #-----------------------------------------------------------------------
        # Axis properties
        #-----------------------------------------------------------------------
        axis.set_ylabel(sqw_labels[proj],fontsize=Correlation_Plots.font_size)
        axis.set_xticks([])
        axis.tick_params(axis='x', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        axis.tick_params(axis='y', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        for ax in ['top','bottom','left','right']:
            axis.spines[ax].set_linewidth(3)
        #-----------------------------------------------------------------------
        # Plotting the S(q,w)
        #-----------------------------------------------------------------------
        axis.imshow(self.sqw_convolve(sqw_data[proj]),origin='lower',aspect='auto',interpolation='quadric',\
        cmap=Correlation_Plots.cmap[col_indx],extent=ax_limits)
        axis.set_xlim(ax_limits[0],ax_limits[1])
        axis.set_ylim(ax_limits[2],ax_limits[3])
        if self.grid and len(q_idx)>0 and len(q_labels)>0:
            axis.set_xticks(q_idx,q_labels)
            axis.xaxis.grid(visible=self.grid,which='major',color='k')
            axis.set_xlabel('',fontsize=Correlation_Plots.font_size)
        else:
            axis.xaxis.grid(visible=self.grid)
            axis.set_xlabel('q',fontsize=Correlation_Plots.font_size)
        return
    ############################################################################
    # @brief Function to plot the spin-spin correlation function obtained from the
    # sqw.simid.out file as well as the AMS obtained from the ams.simid.out file
    # containing the linear spin-wave theory dispersion relation.
    #
    # @author Jonathan Chico
    ############################################################################
    def AMS_Sqw_Plot(self,axis,sqw_data,proj,sqw_labels,ams_data_x,ams_data_y,hf_scale,\
    col_indx,ax_limits,q_labels,q_idx):
        import numpy as np
        axis.cla()
        # Energy maximim 
        self.w_max=ax_limits[3]
        self.w_min=self.w_max/sqw_data[proj].shape[0]
        # Setting the Gaussian smearing to minimum relevant value (0 is not good..)
        self.sigma_w=np.maximum(self.sigma_w,self.w_max/sqw_data[proj].shape[0])
        #-----------------------------------------------------------------------
        # Plotting the S(q,w) and AMS
        #-----------------------------------------------------------------------
        axis.imshow(self.sqw_convolve(sqw_data[proj]),origin='lower',extent=ax_limits,aspect='auto',\
        interpolation='quadric',cmap=Correlation_Plots.cmap[col_indx])
        axis.set_xlim(ax_limits[0],ax_limits[1])
        axis.set_ylim(ax_limits[2],ax_limits[3])
        for ii in range(len(ams_data_y)):
            ams_data_x[ii]=np.arange(ax_limits[0],ax_limits[1]+1,1)
            if ii==0:
                axis.plot(ams_data_x[ii],ams_data_y[ii], c='b',lw=3,label='AMS')
            else:
                axis.plot(ams_data_x[ii],ams_data_y[ii], c='b',lw=3)
        #-----------------------------------------------------------------------
        # Axis properties
        #-----------------------------------------------------------------------
        axis.set_xticks([])
        axis.set_xlabel('q',fontsize=Correlation_Plots.font_size)
        axis.set_ylabel(sqw_labels[proj],fontsize=Correlation_Plots.font_size)
        axis.tick_params(axis='x', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        axis.tick_params(axis='y', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        for ax in ['top','bottom','left','right']:
            axis.spines[ax].set_linewidth(3)
        axis.legend(fontsize=Correlation_Plots.font_size)

        if self.grid and len(q_idx)>0 and len(q_labels)>0:
            axis.set_xticks(q_idx,q_labels)
            axis.xaxis.grid(visible=self.grid,which='major',color='k')
            axis.set_xlabel('',fontsize=Correlation_Plots.font_size)
        else:
            axis.xaxis.grid(visible=self.grid)
        return
