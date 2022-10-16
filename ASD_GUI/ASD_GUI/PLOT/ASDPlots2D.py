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
        return
    ############################################################################
    # @brief Function to plot line and scatter plots.
    # @details This function takes an array of arrays and plot them sequentially
    # the scatter points are located on top of solid lines, with colors been given
    # by a Paired colormap whose size is determined by the size of the array.
    # @author Jonathan Chico
    ############################################################################
    def LinePlot(self,axis,data_x,data_y,labels,ax_label):
        from matplotlib import cm as cm
        import numpy as np
        axis.cla()
        # AB -> add figure properties
        colors=cm.Paired(np.linspace(0,1,len(data_x)+2))
        print('AB  ->',Abstract2DPlot.markersize,self.markersize)
        for ii in range(0,len(data_x)):
            axis.plot(data_x[ii],data_y[ii],lw=self.linewidth,c=colors[ii],label=labels[ii],zorder=-1,markersize=self.markersize,marker='o')
            #axis.scatter(data_x[ii],data_y[ii],color=colors[ii],alpha=0.75, s=150,lw=1.00, edgecolor='black')
        axis.set_xlabel(ax_label[0],fontsize=Abstract2DPlot.font_size)
        axis.set_ylabel(ax_label[1],fontsize=Abstract2DPlot.font_size)
        axis.tick_params(axis='x', colors='black',labelsize=Abstract2DPlot.font_size,width=2)
        axis.tick_params(axis='y', colors='black',labelsize=Abstract2DPlot.font_size,width=2)
        axis.legend(fontsize=Abstract2DPlot.font_size)
        for ax in ['top','bottom','left','right']:
            axis.spines[ax].set_linewidth(3)
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
        return
    ############################################################################
    # @brief Function to plot the spin-spin correlation function obtained from the
    # sqw.simid.out file.
    #
    # @author Jonathan Chico
    ############################################################################
    def Sqw_Plot(self,axis,sqw_data,proj,sqw_labels,col_indx,ax_limits):
        axis.cla()
        #-----------------------------------------------------------------------
        # Axis properties
        #-----------------------------------------------------------------------
        axis.set_xlabel('q',fontsize=Correlation_Plots.font_size)
        axis.set_ylabel(sqw_labels[proj],fontsize=Correlation_Plots.font_size)
        axis.set_xticks([])
        axis.tick_params(axis='x', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        axis.tick_params(axis='y', colors='black',labelsize=Correlation_Plots.font_size,width=2)
        for ax in ['top','bottom','left','right']:
            axis.spines[ax].set_linewidth(3)
        #-----------------------------------------------------------------------
        # Plotting the S(q,w)
        #-----------------------------------------------------------------------
        axis.imshow(sqw_data[proj],origin='lower',aspect='auto',interpolation='lanczos',\
        cmap=Correlation_Plots.cmap[col_indx],extent=ax_limits)
        axis.set_xlim(ax_limits[0],ax_limits[1])
        axis.set_ylim(ax_limits[2],ax_limits[3])
        return
    ############################################################################
    # @brief Function to plot the spin-spin correlation function obtained from the
    # sqw.simid.out file as well as the AMS obtained from the ams.simid.out file
    # containing the linear spin-wave theory dispersion relation.
    #
    # @author Jonathan Chico
    ############################################################################
    def AMS_Sqw_Plot(self,axis,sqw_data,proj,sqw_labels,ams_data_x,ams_data_y,hf_scale,\
    col_indx,ax_limits):
        import numpy as np
        axis.cla()
        #-----------------------------------------------------------------------
        # Plotting the S(q,w) and AMS
        #-----------------------------------------------------------------------
        axis.imshow(sqw_data[proj],origin='lower',extent=ax_limits,aspect='auto',\
        interpolation='quadric',cmap=Correlation_Plots.cmap[col_indx])
        axis.set_xlim(ax_limits[0],ax_limits[1])
        axis.set_ylim(ax_limits[2],ax_limits[3])
        for ii in range(len(ams_data_y)):
            ams_data_x[ii]=np.arange(ax_limits[0]+1,ax_limits[1],1)
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
        return
