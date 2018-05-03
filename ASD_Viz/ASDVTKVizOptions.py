#!/usr/bin/env vtkpython
################################################################################
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper class for the visualization options
################################################################################
import vtk
import ASDMomVTKActors
import ASDNeighVTKActors

class ASDVizOptions():

	MomActors=ASDMomVTKActors.ASDMomActors()
	NeighActors=ASDNeighVTKActors.ASDNeighActors()

	############################################################################
	# Set the color map to be given by the diverging Coolwarm scheme by Kenneth Moreland
	############################################################################
	def set_Coolwarm_lut(self,check):

	    if check:
	        if ASDVizOptions.MomActors.glob_flag_2D:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
	            ASDVizOptions.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
	            ASDVizOptions.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.MagDensMap.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
	        else:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
	            ASDVizOptions.transfer_func.AddRGBPoint(-1, 0.230, 0.299, 0.754)
	            ASDVizOptions.transfer_func.AddRGBPoint( 1, 0.706, 0.016, 0.150)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)

	############################################################################
	# Set the color to be given by the black body function
	############################################################################
	def set_BlackBody_lut(check):

		if check:
			if ASDVizOptions.MomActors.glob_flag_2D:
				ASDVizOptions.lut = vtk.vtkLookupTable()
				num_colors = 256
				ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
				ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
				ASDVizOptions.transfer_func.SetColorSpaceToRGB();
				ASDVizOptions.transfer_func.AddRGBPoint(0.0, 0.0, 0.0, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint(0.4, 0.9, 0.0, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint(0.8, 0.9, 0.9, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint(1.0, 1.0, 1.0, 1.0);
				for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
					cc = ASDVizOptions.transfer_func.GetColor(ss)
					ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
				ASDVizOptions.lut.Build()
				ASDVizOptions.MomActors.MagDensMap.SetLookupTable(ASDVizOptions.lut)
				ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.lut)
				ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
				ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
			else:
				ASDVizOptions.lut = vtk.vtkLookupTable()
				num_colors = 256
				ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
				ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
				ASDVizOptions.transfer_func.SetColorSpaceToRGB();
				ASDVizOptions.transfer_func.AddRGBPoint(-1.0, 0.0, 0.0, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint(-0.5, 0.9, 0.0, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint( 0.5, 0.9, 0.9, 0.0);
				ASDVizOptions.transfer_func.AddRGBPoint( 1.0, 1.0, 1.0, 1.0);
				for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
					cc = ASDVizOptions.transfer_func.GetColor(ss)
					ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
				ASDVizOptions.lut.Build()
				ASDVizOptions.MomActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
				ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
				ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.transfer_func)
				ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)
	############################################################################
	# Set the color map to be given by the diverging RdGy
	############################################################################
	def set_RdGy_lut(self,check):

	    if check:
	        if ASDVizOptions.MomActors.glob_flag_2D:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
	            ASDVizOptions.transfer_func.AddRGBPoint(0.0, 0.79216, 0.00000, 0.12549)
	            ASDVizOptions.transfer_func.AddRGBPoint(0.5, 1.00000, 1.00000, 1.00000)
	            ASDVizOptions.transfer_func.AddRGBPoint(1.0, 0.25098, 0.25098, 0.25098)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.MagDensMap.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
	        else:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
	            ASDVizOptions.transfer_func.AddRGBPoint(-1.0, 0.79216, 0.00000, 0.12549)
	            ASDVizOptions.transfer_func.AddRGBPoint( 0.0, 1.00000, 1.00000, 1.00000)
	            ASDVizOptions.transfer_func.AddRGBPoint( 1.0, 0.25098, 0.25098, 0.25098)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)

	############################################################################
	# Set the color map to be given by the diverging spectral clor map
	############################################################################
	def set_Spectral_lut(self,check):

	    if check:
	        if ASDVizOptions.MomActors.glob_flag_2D:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToRGB()
	            ASDVizOptions.transfer_func.AddRGBPoint(0.00, 0.61961, 0.00392, 0.25882)
	            ASDVizOptions.transfer_func.AddRGBPoint(0.25, 0.95686, 0.42745, 0.26275)
	            ASDVizOptions.transfer_func.AddRGBPoint(0.50, 1.00000, 1.00000, 0.74902)
	            ASDVizOptions.transfer_func.AddRGBPoint(0.75, 0.40000, 0.76078, 0.64706)
	            ASDVizOptions.transfer_func.AddRGBPoint(1.00, 0.36863, 0.30980, 0.63529)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.MagDensMap.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
	        else:
	            ASDVizOptions.lut = vtk.vtkLookupTable()
	            num_colors = 256
	            ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
	            ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
	            ASDVizOptions.transfer_func.SetColorSpaceToRGB()
	            ASDVizOptions.transfer_func.AddRGBPoint(-1.00, 0.61961, 0.00392, 0.25882)
	            ASDVizOptions.transfer_func.AddRGBPoint(-0.50, 0.95686, 0.42745, 0.26275)
	            ASDVizOptions.transfer_func.AddRGBPoint( 0.00, 1.00000, 1.00000, 0.74902)
	            ASDVizOptions.transfer_func.AddRGBPoint( 0.50, 0.40000, 0.76078, 0.64706)
	            ASDVizOptions.transfer_func.AddRGBPoint( 1.00, 0.36863, 0.30980, 0.63529)
	            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	                cc = ASDVizOptions.transfer_func.GetColor(ss)
	                ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	            ASDVizOptions.lut.Build()
	            ASDVizOptions.MomActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.transfer_func)
	            ASDVizOptions.MomActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)

	############################################################################
	# Toggle option for the axes
	############################################################################
	def toggle_Axes(self,check):
		if check:
			ASDVizOptions.MomActors.OrientMarker.SetEnabled(1)
		else:
			ASDVizOptions.MomActors.OrientMarker.SetEnabled(0)

	############################################################################
	# Toggle option for the scalar bar
	############################################################################
	def toggle_ScalarBar(self,check):
	    if check:
	        ASDVizOptions.MomActors.scalar_bar_widget.SetEnabled(1)
	    else:
	        ASDVizOptions.MomActors.scalar_bar_widget.SetEnabled(0)

	############################################################################
	# Toggle options for the contours
	############################################################################
	def toggle_contours(self,check):
	    if check:
	        ASDVizOptions.MomActors.contActor.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.contActor.VisibilityOff()

	############################################################################
	# Toggle the directions arrows
	############################################################################
	def toggle_directions(self,check):
	    if check:
	        ASDVizOptions.MomActors.vector.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.vector.VisibilityOff()

	############################################################################
	# Toggle the directions arrows
	############################################################################
	def toggle_spins(self,check):
	    if check:
	        ASDVizOptions.MomActors.Spins.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.Spins.VisibilityOff()

	############################################################################
	# Toggle the magnetization density
	############################################################################
	def toggle_density(self,check):
	    if check:
	        ASDVizOptions.MomActors.MagDensActor.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.MagDensActor.VisibilityOff()

	############################################################################
	# Toggle the visualization of the embeded cluster
	############################################################################
	def toggle_cluster(self,check):
	    if check:
	        ASDVizOptions.MomActors.atom.VisibilityOn()
	        ASDVizOptions.MomActors.atom_imp.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.atom.VisibilityOff()
	        ASDVizOptions.MomActors.atom_imp.VisibilityOff()

	############################################################################
	# Toggle the KMC particle visualization
	############################################################################
	def toggle_KMC(self,check):
	    if check:
	        ASDVizOptions.MomActors.KMC_part_actor.VisibilityOn()
	    else:
	        ASDVizOptions.MomActors.KMC_part_actor.VisibilityOff()

	############################################################################
	# Toggle the plane clipper
	############################################################################
	def toggle_clipper(self,check):
	    if check:
	        ASDVizOptions.MomActors.clipperActor.VisibilityOn()
	        ASDVizOptions.MomActors.MagDensActor.VisibilityOff()
	    else:
	        ASDVizOptions.MomActors.clipperActor.VisibilityOff()
	        ASDVizOptions.MomActors.MagDensActor.VisibilityOn()

	############################################################################
	# Set the color of the magnetization density along the x axis projection
	############################################################################
	def set_color_x (self,check):
	    if check:
	        ASDVizOptions.MomActors.src.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_x)

	############################################################################
	# Set the color of the magnetization density along the y axis projection
	############################################################################
	def set_color_y (self,check):
	    if check:
	        ASDVizOptions.MomActors.src.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_y)

	############################################################################
	# Set the color of the magnetization density along the z axis projection
	############################################################################
	def set_color_z (self,check):
	    if check:
	        ASDVizOptions.MomActors.src.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_z)

	############################################################################
	# Set the color of the magnetization density along the x axis projection for the spins
	############################################################################
	def set_spins_color_x (self,check):
	    if check:
	        ASDVizOptions.MomActors.src_spins.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_x)

	############################################################################
	# Set the color of the magnetization density along the y axis projection
	############################################################################
	def set_spins_color_y (self,check):
	    if check:
	        ASDVizOptions.MomActors.src_spins.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_y)

	############################################################################
	# Set the color of the magnetization density along the z axis projection
	############################################################################
	def set_spins_color_z (self,check):
	    if check:
	        ASDVizOptions.MomActors.src_spins.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color_z)

	############################################################################
	# Set the size of the spins via the slider
	############################################################################
	def ChangeSpinsSize (self,value):
	    ASDVizOptions.MomActors.SpinMapper.SetScaleFactor(0.50*value/10)

    ############################################################################
    # Toggle the atoms for the neighbour map
    ############################################################################
	def toggle_NAtoms(self,check):
		if check:
			ASDVizOptions.NeighActors.AtomsActor.VisibilityOn()
		else:
			ASDVizOptions.NeighActors.AtomsActor.VisibilityOff()

    ############################################################################
    # Toggle the neighbour cloud for the neighbour map
    ############################################################################
	def toggle_Neigh(self,check):
		if check:
			ASDVizOptions.NeighActors.NeighActor.VisibilityOn()
		else:
			ASDVizOptions.NeighActors.NeighActor.VisibilityOff()
