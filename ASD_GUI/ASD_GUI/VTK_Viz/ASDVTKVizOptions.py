"""@package ASDVTKVizOptions
Contains a class with a set of functions dealing with the visualization options
in the VTK visualization mode. Dealing with visibility of objects, size of glyphs,
types of glyphs, colormap used etc.

Author
----------
Jonathan Chico
"""
from ASD_GUI.VTK_Viz import ASDVTKMomActors 
from ASD_GUI.VTK_Viz import ASDVTKGenActors 
from ASD_GUI.VTK_Viz import ASDVTKEneActors 
from ASD_GUI.VTK_Viz import ASDVTKNeighActors 

################################################################################
# @brief Class containing the majority of the actions to update the visualizer
# @details Class containing the majority of the actions to update the visualizer.
# It handles the taking of snapshots, the toggling of projections, actor visibility
# and colormaps. It also controls many of the camera options so that they can be
# updated during runtime, allowing for finer control of the visualization.
# @author Jonathan Chico
################################################################################
class ASDVizOptions():
    from vtkmodules.vtkCommonCore import vtkLookupTable 
    GenActors=ASDVTKGenActors.ASDGenActors()
    EneActors=ASDVTKEneActors.ASDEneActors()
    MomActors=ASDVTKMomActors.ASDMomActors()
    NeighActors=ASDVTKNeighActors.ASDNeighActors()

    lut = vtkLookupTable()

    ############################################################################
    # @ brief A function that takes a renderwindow and saves its contents to a .png file
    # @author Anders Bergman
    ############################################################################
    def Screenshot(self,renWin,number_of_screenshots,png_mode,pov_mode):
        """Function to take the rendering window and save it to file, either in .png
        or in .pov format.

        Args:
            renWin: current rendering window.
            number_of_screenshots: current number of the screenshot that is being saved.
            png_mode: (logical) variable indicated if the scene should be stored as a png
            pov_mode: (logical) variable indicated if the scene should be stored as a pov

        Author
        ----------
        Anders Bergman
        """
        from vtk import vtkWindowToImageFilter,vtkPOVExporter,vtkPNGWriter

        win2im=vtkWindowToImageFilter()
        win2im.SetInput(renWin)
        win2im.Update()
        win2im.SetInputBufferTypeToRGBA()
        win2im.ReadFrontBufferOff()
        #-----------------------------------------------------------------------
        # Save snapshot as a '.pov'
        #-----------------------------------------------------------------------
        if pov_mode:
            povexp=vtkPOVExporter()
            povexp.SetInput(renWin)
            renWin.Render()
            povexp.SetFileName('snap%.5d.pov' %number_of_screenshots)
            povexp.Write()
        #-----------------------------------------------------------------------
        # Save snapshot as a '.png'
        #-----------------------------------------------------------------------
        if png_mode:
            toPNG=vtkPNGWriter()
            toPNG.SetFileName('snap%.5d.png' %number_of_screenshots)
            toPNG.SetInputConnection(win2im.GetOutputPort())
            toPNG.Write()
        return;
    ############################################################################
    # @brief Toggle the parallel projection for the camera
    # @author Jonathan Chico
    ############################################################################
    def toggle_projections(self,renWin,window,ren,checked):
        if checked:
            ren.GetActiveCamera().ParallelProjectionOn()
            window.ParallelScaleLineEdit.setText(str(window.ParallelScaleSlider.value()))
            ren.GetActiveCamera().SetParallelScale(float(window.ParallelScaleLineEdit.text()))
            renWin.Render()
        else:
            ren.GetActiveCamera().ParallelProjectionOff()
            renWin.Render()
        return
    ############################################################################
    # @brief Change the parameters of the parallel projection
    # @author Jonathan Chico
    ############################################################################
    def ChangeParallelProj(self,ren,renWin,line,slider,MainWindow):
        if line:
            MainWindow.ParallelScaleSlider.setValue(float(MainWindow.ParallelScaleLineEdit.text()))
            ren.GetActiveCamera().SetParallelScale(float(MainWindow.ParallelScaleLineEdit.text()))
            renWin.Render()
        if slider:
            ren.GetActiveCamera().SetParallelScale(MainWindow.ParallelScaleSlider.value())
            MainWindow.ParallelScaleLineEdit.setText(str(MainWindow.ParallelScaleSlider.value()))
            renWin.Render()
        return
    ############################################################################
    # @brief Function to reset the camera to the initial positions
    # @author Jonathan Chico
    ############################################################################
    def reset_camera(self,ren,renWin,current_Actors):
        """Function to reset the camera to the initial position.

        Args:
            ren: current VTK renderer.
            renWin: current VTK rendering window.
            current_Actors: current actors which are being visualized.

        Author
        ----------
        Jonathan Chico
        """
        # Defining the camera directions
        ren.GetActiveCamera().SetFocalPoint(current_Actors.xmid,current_Actors.ymid,current_Actors.zmid)
        ren.GetActiveCamera().SetPosition(current_Actors.xmid,current_Actors.ymid,current_Actors.height)
        ren.GetActiveCamera().Azimuth(0)
        ren.GetActiveCamera().Elevation(0)
        ren.GetActiveCamera().Yaw(0)
        ren.GetActiveCamera().Roll(0)
        ren.GetActiveCamera().Pitch(0)
        ren.GetActiveCamera().SetViewUp(0,1,0)
        renWin.Render()
        return
    ############################################################################
    # @brief Update the camera values for those defined by the user in the GUI.
    # @author Jonathan Chico
    ############################################################################
    def Update_Camera(self,Window,ren,renWin):
        """Function to update the value of the camera for the current visualization from
        the values entered by the user in the GUI.

        Args:
            Window: QMainWindow where the visualizations are being carried out.
            ren: current VTK renderer.
            renWin: current VTK rendering window.

        Author
        ----------
        Jonathan Chico
        """
        camera_focal=[0]*3
        camera_pos=[0]*3
        camera_pos[0]=float(Window.CamPosX.text())
        camera_pos[1]=float(Window.CamPosY.text())
        camera_pos[2]=float(Window.CamPosZ.text())
        camera_focal[0]=float(Window.FocalPosX.text())
        camera_focal[1]=float(Window.FocalPosY.text())
        camera_focal[2]=float(Window.FocalPosZ.text())
        ren.GetActiveCamera().SetFocalPoint(camera_focal)
        ren.GetActiveCamera().SetPosition(camera_pos)
        ren.GetActiveCamera().Elevation(float(Window.CamElevationLineEdit.text()))
        ren.GetActiveCamera().Azimuth(float(Window.CamAzimuthLineEdit.text()))
        ren.GetActiveCamera().Pitch(float(Window.CamPitchLineEdit.text()))
        ren.GetActiveCamera().Roll(float(Window.CamRollLineEdit.text()))
        ren.GetActiveCamera().Yaw(float(Window.CamYawLineEdit.text()))
        renWin.Render()
        return
    ############################################################################
    # Update the dock window information when the camera is setup
    ############################################################################
    def update_dock_info(self,current_Actors,Window):
        Window.FocalPosX.setText(str(current_Actors.camera_focal[0]))
        Window.FocalPosY.setText(str(current_Actors.camera_focal[1]))
        Window.FocalPosZ.setText(str(current_Actors.camera_focal[2]))
        Window.CamPosX.setText(str(current_Actors.camera_pos[0]))
        Window.CamPosY.setText(str(current_Actors.camera_pos[1]))
        Window.CamPosZ.setText(str(current_Actors.camera_pos[2]))
        Window.CamYawLineEdit.setText(str(current_Actors.camera_yaw))
        Window.CamRollLineEdit.setText(str(current_Actors.camera_roll))
        Window.CamPitchLineEdit.setText(str(current_Actors.camera_pitch))
        Window.CamAzimuthLineEdit.setText(str(current_Actors.camera_azimuth))
        Window.CamElevationLineEdit.setText(str(current_Actors.camera_elevation))
        return
    ############################################################################
    # Set the camera view up to be defined to the (1,0,0)
    ############################################################################
    def set_Camera_viewUp(self,ren,renWin,dir):
        ren.GetActiveCamera().SetViewUp(dir)
        renWin.Render()
        return
    ############################################################################
    # Toggle option for the axes
    ############################################################################
    def toggle_Axes(self,check):
        if check:
            ASDVizOptions.GenActors.OrientMarker.SetEnabled(1)
        else:
            ASDVizOptions.GenActors.OrientMarker.SetEnabled(0)
        return
    ############################################################################
    # Toggle option for the scalar bar
    ############################################################################
    def toggle_ScalarBar(self,check):
        if check:
            ASDVizOptions.GenActors.scalar_bar_widget.SetEnabled(1)
        else:
            ASDVizOptions.GenActors.scalar_bar_widget.SetEnabled(0)
        return
    ############################################################################
    # Toggle options for the contours
    ############################################################################
    def toggle_contours(self,check):
        if check:
            ASDVizOptions.MomActors.contActor.VisibilityOn()
        else:
            ASDVizOptions.MomActors.contActor.VisibilityOff()
        return
    ############################################################################
    # Toggle the directions arrows
    ############################################################################
    def toggle_directions(self,check):
        if check:
            ASDVizOptions.MomActors.vector.VisibilityOn()
        else:
            ASDVizOptions.MomActors.vector.VisibilityOff()
        return
    ############################################################################
    # Toggle the directions arrows
    ############################################################################
    def toggle_spins(self,check):
        if check:
            ASDVizOptions.MomActors.Spins.VisibilityOn()
        else:
            ASDVizOptions.MomActors.Spins.VisibilityOff()
        return
    ############################################################################
    # Toggle the atomic spheres
    ############################################################################
    def toggle_atoms(self,check):
        if check:
            ASDVizOptions.MomActors.Atoms.VisibilityOn()
        else:
            ASDVizOptions.MomActors.Atoms.VisibilityOff()
        return
    ############################################################################
    # Toggle the magnetization density
    ############################################################################
    def toggle_density(self,check):
        if check:
            ASDVizOptions.MomActors.MagDensActor.VisibilityOn()
        else:
            ASDVizOptions.MomActors.MagDensActor.VisibilityOff()
        return
    ############################################################################
    # Toggle the visualization of the embedded cluster
    ############################################################################
    def toggle_cluster(self,check):
        if check:
            ASDVizOptions.GenActors.atom.VisibilityOn()
            ASDVizOptions.GenActors.atom_imp.VisibilityOn()
        else:
            ASDVizOptions.GenActors.atom.VisibilityOff()
            ASDVizOptions.GenActors.atom_imp.VisibilityOff()
        return
    ############################################################################
    # Toggle the KMC particle visualization
    ############################################################################
    def toggle_KMC(self,check):
        if check:
            ASDVizOptions.MomActors.KMC_part_actor.VisibilityOn()
        else:
            ASDVizOptions.MomActors.KMC_part_actor.VisibilityOff()
        return
    ############################################################################
    # Toggle the plane clipper
    ############################################################################
    def toggle_clipper(self,check,current_Actors,dir,window,origin,min,max,renWin):
        if check:
            ASDVizOptions.GenActors.clipperActor.VisibilityOn()
            current_Actors.VisibilityOff()
            self.set_clipp_plane(dir,window,origin,min,max,renWin)
        else:
            ASDVizOptions.GenActors.clipperActor.VisibilityOff()
            current_Actors.VisibilityOn()
            renWin.Render()
        return
    ############################################################################
    # Toggle the time label
    ############################################################################
    def toggle_time_label(self,check):
        if check:
            ASDVizOptions.GenActors.time_label_widget.On()
        else:
            ASDVizOptions.GenActors.time_label_widget.Off()
        return
    ############################################################################
    # Set the clipping plane such that the normal plane is 'origin'
    ############################################################################
    def set_clipp_plane(self,dir,window,origin,min,max,renWin):
        ASDVizOptions.GenActors.plane.SetOrigin(origin)
        ASDVizOptions.GenActors.plane.SetNormal(dir)
        window.ClippingPlaneSlider.setMinimum(min)
        window.ClippingPlaneSlider.setMaximum(max)
        window.ClippingPlaneSlider.setValue(min)
        renWin.Render()
        return
    ############################################################################
    # Set the position of the clipping plane via the slider
    ############################################################################
    def ClippingUpdate(self,origin,window,renWin):
        ASDVizOptions.GenActors.plane.SetOrigin(origin)
        window.ClipPlaneLabel.setText('Clip. Plane Pos.={:.1f},{:.1f},{:.1f}'\
        .format(float(origin[0]),float(origin[1]),float(origin[2])))
        renWin.Render()
        return
    ############################################################################
    # Set the color of the magnetization density along the x axis projection
    ############################################################################
    def set_projection(self,type,axis):
        if type=='density':
            ASDVizOptions.MomActors.src.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color[axis])
        elif type=='spins':
            ASDVizOptions.MomActors.src_spins.GetPointData().SetScalars(ASDVizOptions.MomActors.glob_color[axis])
        return
    ############################################################################
    # Set the size of the spins via the slider
    ############################################################################
    def ChangeSpinsSize(self,value):
        ASDVizOptions.MomActors.SpinMapper.SetScaleFactor(0.50*value/10)
        return

    ############################################################################
    # Set the interpolation of the spin glyphs
    ############################################################################
    def ChangeSpinShade(self,renWin,keyword):
        if keyword=='Flat':
            if hasattr(ASDVizOptions.MomActors,'Spins'):
                ASDVizOptions.MomActors.Spins.GetProperty().SetInterpolationToFlat()
        elif keyword=='Gouraud':
            if hasattr(ASDVizOptions.MomActors,'Spins'):
                ASDVizOptions.MomActors.Spins.GetProperty().SetInterpolationToGouraud()
            if hasattr(ASDVizOptions.MomActors,'Atoms'):
                ASDVizOptions.MomActors.Atoms.GetProperty().SetInterpolationToGouraud()
        elif keyword=='PBR':
            if hasattr(ASDVizOptions.MomActors,'Spins'):
                ASDVizOptions.MomActors.Spins.GetProperty().SetInterpolationToPBR()
                ASDVizOptions.MomActors.Spins.GetProperty().SetMetallic(0.5)
            if hasattr(ASDVizOptions.MomActors,'Atoms'):
                ASDVizOptions.MomActors.Atoms.GetProperty().SetInterpolationToPBR()
                ASDVizOptions.MomActors.Atoms.GetProperty().SetMetallic(0.5)
        elif keyword=='Phong':
            if hasattr(ASDVizOptions.MomActors,'Spins'):
                ASDVizOptions.MomActors.Spins.GetProperty().SetInterpolationToPhong()
            if hasattr(ASDVizOptions.MomActors,'Atoms'):
                ASDVizOptions.MomActors.Atoms.GetProperty().SetInterpolationToPhong()

        renWin.Render()
        return

    ############################################################################
    # Set the ambient scattering of the spin glyphs
    ############################################################################
    def RenAmbientUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetAmbient(float(value*0.02))
        renWin.Render()
        return
    ############################################################################
    # Set the diffuse scattering of the spin glyphs
    ############################################################################
    def RenDiffuseUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetDiffuse(float(value*0.01))
        renWin.Render()
        return
    ############################################################################
    # Set the specular scattering of the spin glyphs
    ############################################################################
    def RenSpecularUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetSpecular(float(value*0.01))
        renWin.Render()
        return

    ############################################################################
    # Set the specular scattering of the spin glyphs
    ############################################################################
    def RenSpecularPowerUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetSpecularPower(float(value))
        renWin.Render()
        return

    ############################################################################
    # Set the PBR emission value of the spin glyphs
    ############################################################################
    def PBREmissionUpdate(self,value,ren, renWin):
        emvec = [float(value*0.01), float(value*0.01), float(value*0.01)]
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetEmissiveFactor(emvec)
        if hasattr(ASDVizOptions.MomActors, 'Atoms'):
            ASDVizOptions.MomActors.Atoms.GetProperty().SetEmissiveFactor(emvec)
        renWin.Render()
        return
    ############################################################################
    # Set the PBR occlusion value of the spin glyphs
    ############################################################################
    def PBROcclusionUpdate(self,value,ren, renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetOcclusionStrength(float(value*0.01))
        if hasattr(ASDVizOptions.MomActors, 'Atoms'):
            ASDVizOptions.MomActors.Atoms.GetProperty().SetOcclusionStrength(float(value*0.01))
        renWin.Render()
        return
    ############################################################################
    # Set the PBR roughness value of the spin glyphs
    ############################################################################
    def PBRRoughnessUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetRoughness(float(value*0.01))
        if hasattr(ASDVizOptions.MomActors, 'Atoms'):
            ASDVizOptions.MomActors.Atoms.GetProperty().SetRoughness(float(value*0.01))
        renWin.Render()
        return
    ############################################################################
    # Set the PBR metallic value of the spin glyphs
    ############################################################################
    def PBRMetallicUpdate(self,value,renWin):
        if hasattr(ASDVizOptions.MomActors, 'Spins'):
            ASDVizOptions.MomActors.Spins.GetProperty().SetMetallic(float(value*0.01))
        if hasattr(ASDVizOptions.MomActors, 'Atoms'):
            ASDVizOptions.MomActors.Atoms.GetProperty().SetMetallic(float(value*0.01))
        renWin.Render()
        return

    ############################################################################
    # Set the size of the atoms via the slider
    ############################################################################
    def ChangeAtomsSize(self,value):
        ASDVizOptions.MomActors.AtomMapper.SetScaleFactor(1.00*value/10.0)
        return
    ############################################################################
    # Set the size of the atoms via the slider
    ############################################################################
    def ChangeAtomsOpaq(self,value):
        ASDVizOptions.MomActors.Atoms.GetProperty().SetOpacity(value*0.01)
        return
    ############################################################################
    # Toggle the atoms for the neighbour map
    ############################################################################
    def toggle_NAtoms(self,check):
        if check:
            ASDVizOptions.NeighActors.AtomsActor.VisibilityOn()
        else:
            ASDVizOptions.NeighActors.AtomsActor.VisibilityOff()
        return
    ############################################################################
    # Toggle the neighbour cloud for the neighbour map
    ############################################################################
    def toggle_Neigh(self,check):
        if check:
            ASDVizOptions.NeighActors.NeighActor.VisibilityOn()
        else:
            ASDVizOptions.NeighActors.NeighActor.VisibilityOff()
        return
    ############################################################################
    # Set the opacity of the neighbour spheres
    ############################################################################
    def NeighOpacityUpdate(self,value):
        ASDVizOptions.NeighActors.NeighActor.GetProperty().SetOpacity(value*0.1)
        return
    ############################################################################
    # Set the opacity of the atom spheres
    ############################################################################
    def AtomOpacityUpdate(self,value):
        ASDVizOptions.NeighActors.AtomsActor.GetProperty().SetOpacity(value*0.1)
        return
    ############################################################################
    ## @brief Function to find the needed file names for the HDR file
    ############################################################################
    def getHDRIFileName(self,window):
        from PyQt6 import QtWidgets

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
        hdrifile = dlg.getOpenFileName(caption="Open HDR file",
                        directory= '.', filter="HDR images (*.pic *.hdr)")[0]
        return hdrifile
    ############################################################################
    ## @brief Function to find the needed file names for the texture files
    ############################################################################
    def getTextureFileName(self,window):
        from PyQt6 import QtWidgets

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
        texturefile = dlg.getOpenFileName(caption="Open texture file", 
                                directory= '.', filter="Images (*.png)")[0]
        return texturefile
    ############################################################################
    # Toggle surface texture 
    ############################################################################
    def toggle_Texture(self,check, ren, renWin, texfile):
        import vtk

        if check:

            albedo_reader = vtk.vtkPNGReader()
            albedo_reader.SetFileName(texfile)

            albedo = vtk.vtkTexture()
            albedo.SetInputConnection(albedo_reader.GetOutputPort())
            albedo.UseSRGBColorSpaceOn()
            albedo.InterpolateOn()
            albedo.MipmapOn()

            self.MomActors.Spins.GetProperty().SetBaseColorTexture(albedo)
        else:
            self.MomActors.Spins.GetProperty().RemoveTexture('albedoTex')


        renWin.Render()

        return
    ############################################################################
    # Toggle ORM texture 
    ############################################################################
    def toggle_ORMTexture(self,check, ren, renWin, texfile):
        import vtk

        if check:

            material_reader = vtk.vtkPNGReader()
            material_reader.SetFileName(texfile)

            material = vtk.vtkTexture()
            material.SetInputConnection(material_reader.GetOutputPort())
            material.InterpolateOn()
            material.MipmapOn()

            self.MomActors.Spins.GetProperty().SetORMTexture(material)
        else:
            self.MomActors.Spins.GetProperty().RemoveTexture('materialTex')

        renWin.Render()

        return
    ############################################################################
    # Toggle anisotropy texture 
    ############################################################################
    def toggle_ATexture(self,check, ren, renWin, texfile):
        import vtk

        if check:

            anisotropy_reader = vtk.vtkPNGReader()
            anisotropy_reader.SetFileName(texfile)

            anisotropy = vtk.vtkTexture()
            anisotropy.SetInputConnection(anisotropy_reader.GetOutputPort())
            anisotropy.InterpolateOn()
            anisotropy.MipmapOn()

            self.MomActors.Spins.GetProperty().SetAnisotropyTexture(anisotropy)

        else:
            self.MomActors.Spins.GetProperty().RemoveTexture('anisotropyTex')

        renWin.Render()

        return
    ############################################################################
    # Toggle normal texture 
    ############################################################################
    def toggle_NTexture(self,check, ren, renWin, texfile):
        import vtk

        if check:

            normal_reader = vtk.vtkPNGReader()
            normal_reader.SetFileName(texfile)

            normal = vtk.vtkTexture()
            normal.InterpolateOn()
            normal.MipmapOn()
            normal.SetInputConnection(normal_reader.GetOutputPort())

            self.MomActors.SpinShader=self.MomActors.Spins.GetShaderProperty()

            self.MomActors.SpinShader.AddVertexShaderReplacement(
                    "//VTK::Normal::Dec", # replace the normal block
                    True, # before the standard replacements
                    "//VTK::Normal::Dec\n" # we still want the default
                    "in vec3 tangentMC;\n"
                    "out vec3 tangentVCVSOutput;\n",
                    False # only do it once
                )
            self.MomActors.SpinShader.AddVertexShaderReplacement(
                    "//VTK::Normal::Impl", # replace the normal block
                    True, # before the standard replacements
                    "//VTK::Normal::Impl\n" # we still want the default
                    "  tangentVCVSOutput = normalMatrix * tangentMC;\n",
                    False # only do it once
                )
            self.MomActors.Spins.GetProperty().SetNormalTexture(normal)

        else:
            self.MomActors.SpinShader.ClearAllVertexShaderReplacements()
            self.MomActors.Spins.GetProperty().RemoveTexture('normalTex')

        renWin.Render()

        return
    ############################################################################
    # Toggle ORM texture 
    ############################################################################
    def toggle_ETexture(self,check, ren, renWin, texfile):
        import vtk

        if check:

            emissive_reader = vtk.vtkPNGReader()
            emissive_reader.SetFileName(texfile)

            emissive = vtk.vtkTexture()
            emissive.SetInputConnection(emissive_reader.GetOutputPort())
            emissive.UseSRGBColorSpaceOn()
            emissive.InterpolateOn()
            emissive.MipmapOn()

            self.MomActors.Spins.GetProperty().SetEmissiveTexture(emissive)


        else:
            self.MomActors.Spins.GetProperty().RemoveTexture('emissiveTex')

        renWin.Render()

        return
    ############################################################################
    # Toggle Skybox on/off
    ############################################################################
    def toggle_SkyBox(self,check, ren, renWin, skyboxfile):
        from vtkmodules.vtkIOImage import vtkHDRReader
        from vtkmodules.vtkRenderingCore import vtkTexture

        if check:
            reader = vtkHDRReader()
            reader.SetFileName(skyboxfile)
            reader.Update()

            texture = vtkTexture()
            texture.InterpolateOn()
            texture.MipmapOn()
            texture.SetColorModeToDirectScalars()
            texture.SetInputConnection(reader.GetOutputPort())

            self.MomActors.SkyBox.SetTexture(texture)
            self.MomActors.SkyBox.SetProjectionToSphere()
            self.MomActors.SkyBox.VisibilityOn()
        else:
            self.MomActors.SkyBox.VisibilityOff()

        renWin.Render()

        return
    ############################################################################
    # Toggle HDRI on/off
    ############################################################################
    def toggle_HDRI(self,check, ren, renWin, hdrifile):
        #import vtk
        from vtkmodules.vtkIOImage import vtkHDRReader
        from vtkmodules.vtkRenderingCore import vtkTexture

        if check:
            reader = vtkHDRReader()
            reader.SetFileName(hdrifile)
            reader.Update()

            texture = vtkTexture()
            texture.InterpolateOn()
            texture.MipmapOn()
            texture.SetColorModeToDirectScalars()
            texture.SetInputConnection(reader.GetOutputPort())

 
            #ren.RemoveAllLights()
            ren.AutomaticLightCreationOff()
            ren.UseImageBasedLightingOn()
            ren.UseSphericalHarmonicsOn()
            ren.SetEnvironmentTexture(texture)

        else:
            ren.UseSphericalHarmonicsOff()
            ren.UseImageBasedLightingOff()
            ren.AutomaticLightCreationOn()

        renWin.Render()
        return
    ############################################################################
    # Toggle SSAO on/off
    ############################################################################
    def toggle_SSAO(self,check, ren):

        if check:
            ren.UseSSAOOn()
            ren.SetSSAOKernelSize(512)
            ren.SetSSAORadius(3.0)
            ren.SetSSAOBias(0.1)
            ren.SSAOBlurOff()

            #self.toggle_HDRI(check=check,ren=ren)

        else:
            ren.UseSSAOOff()

        return
    ############################################################################
    # Toggle automatic focal point determination on/off
    ############################################################################
    def setFocalDisk(self,value, ren, renWin):

        ren.GetActiveCamera().SetFocalDisk(value/200.0)

        renWin.Render()
        
    ############################################################################
    # Toggle automatic focal point determination on/off
    ############################################################################
    def toggle_autoFocus(self,check, renWin):

        if check:
            self.dofPass.AutomaticFocalDistanceOn()
        else:
            self.dofPass.AutomaticFocalDistanceOff()

        renWin.Render()
        return
    ############################################################################
    # Toggle depth of field focus on/off
    ############################################################################
    def toggle_Focus(self,check, ren, renWin):
        import vtk

        if check:
            # create the basic VTK render steps
            self.basicPasses = vtk.vtkRenderStepsPass()

            self.dofPass = vtk.vtkDepthOfFieldPass()
            self.dofPass.AutomaticFocalDistanceOff()
            self.dofPass.SetDelegatePass(self.basicPasses)

            # Tell the renderer to use our render pass pipeline.
            ren.GetActiveCamera().SetFocalDisk(0.5)
            ren.SetPass(self.dofPass)

            renWin.Render()

        else:
            print('DOF render pass can not be disabled.')
            ren.ReleaseGraphicsResources(renWin)


        return
    ############################################################################
    # Toggle FXAA on/off
    ############################################################################
    def toggle_FXAA(self,check, ren, renWin):
        if check:
            ren.UseFXAAOn()
            ren.GetFXAAOptions().SetUseHighQualityEndpoints(True)
            renWin.SetMultiSamples(4)
        else:
            ren.UseFXAAOff()
        return
    ############################################################################
    # Toggle shadows on/off
    ############################################################################
    ### def toggle_Shadows(self,check, ren, renWin):
    ###     from vtkmodules.vtkRenderingOpenGL2 import (
    ###         vtkCameraPass,
    ###         vtkOpaquePass,
    ###         vtkRenderPassCollection,
    ###         vtkSequencePass,
    ###         vtkShadowMapPass
    ###     )
    ###     print('Toggle shadows', check)
    ###     if check:
    ###         seq = vtkSequencePass()
    
    ###         passes = vtkRenderPassCollection()
    
    ###         shadows = vtkShadowMapPass()
    ###         passes.AddItem(shadows.GetShadowMapBakerPass())
    ###         passes.AddItem(shadows)
    
    ###         opaque = vtkOpaquePass()
    ###         passes.AddItem(opaque)
    
    ###         seq.SetPasses(passes)
    
    ###         camera_p = vtkCameraPass()
    ###         camera_p.SetDelegatePass(seq)
    
    ###         # Tell the renderer to use our render pass pipeline.
    ###         ren.SetPass(camera_p)
    ###         renWin.Render()

    ###     return
    ############################################################################
    # Update glyph resolution
    ############################################################################
    def GlyphQualityUpdate(self,value,viz_type,mode,renWin):
        if viz_type=='M':
            try:
                ASDVizOptions.MomActors.spinarrow.SetTipResolution(value)
                ASDVizOptions.MomActors.spinarrow.SetShaftResolution(value)
            except:
                pass
            try:
                ASDVizOptions.MomActors.spinsphere.SetThetaResolution(value)
                ASDVizOptions.MomActors.spinsphere.SetPhiResolution(value)
            except:
                pass
            try:
                ASDVizOptions.MomActors.spincones.SetResolution(value)
            except:
                pass

        if viz_type=='N':
            if mode==1:
                ASDVizOptions.NeighActors.NeighGlyphs.SetThetaResolution(value)
                ASDVizOptions.NeighActors.NeighGlyphs.SetPhiResolution(value)
            if mode==2:
                ASDVizOptions.NeighActors.NeighGlyphs.SetTipResolution(value)
                ASDVizOptions.NeighActors.NeighGlyphs.SetShaftResolution(value)
        if viz_type=='E':
            ASDVizOptions.EneActors.EneAtom.SetThetaResolution(value)
            ASDVizOptions.EneActors.EneAtom.SetPhiResolution(value)

        renWin.Render()
        return
    ############################################################################
    # @brief Change the type of glyphs used to visualize the atomistic spins
    # @details Change the type of glyphs used to visualize the atomistic spins.
    # It allows the user to change the glyphs for the atomistic spins on the fly
    # the user can choose between the following actors:
    #
    #   - Cone
    #   - Arrows
    #   - Spheres
    #   - Cubes
    # @author Jonathan Chico
    ############################################################################
    def ChangeSpinGlyph(self,renWin,keyword):
        """Change the type of glyphs used to visualize the atomistic spins.
        It allows the user to change the glyphs for the atomistic spins on the fly
        the user can choose between the following actors:
            * Cone
            * Arrows
            * Spheres
            * Cubes

        Args:
            renWin: current VTK rendering window.
            keyword: (str) identifier keyword to choose between the different types of glyphs.

        Author
        ----------
        Jonathan Chico
        """

        import vtk

        if keyword=='Cubes':
            try:
                del ASDVizOptions.MomActors.spinarrow
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spinsphere
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincones
            except:
                pass
            ASDVizOptions.MomActors.spincube=vtk.vtkCubeSource()
            ASDVizOptions.MomActors.spincube.SetXLength(1.0)
            ASDVizOptions.MomActors.spincube.SetYLength(1.0)
            ASDVizOptions.MomActors.spincube.SetZLength(1.0)

            ## Calculate TCoords for texturing
            #ASDVizOptions.MomActors.spincubetmap = vtk.vtkTextureMapToSphere()
            #ASDVizOptions.MomActors.spincubetmap.SetInputConnection(ASDVizOptions.MomActors.spincube.GetOutputPort())
            #ASDVizOptions.MomActors.spincubetmap.PreventSeamOn()

            ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spincube.GetOutputPort())
            #ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spincubetmap.GetOutputPort())
            ASDVizOptions.MomActors.SpinMapper.ClampingOn()
            ASDVizOptions.MomActors.SpinMapper.OrientOff()
            renWin.Render()
        if keyword=='Bars':
            try:
                del ASDVizOptions.MomActors.spinarrow
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spinsphere
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincones
            except:
                pass
            ASDVizOptions.MomActors.spincube=vtk.vtkCubeSource()
            ASDVizOptions.MomActors.spincube.SetXLength(2.0)
            ASDVizOptions.MomActors.spincube.SetYLength(0.4)
            ASDVizOptions.MomActors.spincube.SetZLength(0.4)

            ### Calculate TCoords for texturing
            ASDVizOptions.MomActors.spincubetmap = vtk.vtkTextureMapToCylinder()
            ASDVizOptions.MomActors.spincubetmap.SetInputConnection(ASDVizOptions.MomActors.spincube.GetOutputPort())
            #ASDVizOptions.MomActors.spincubetmap.AutomaticCylinderGenerationOff()
            #ASDVizOptions.MomActors.spincubetmap.SetPoint1([ 1.0,0.0,0.0])
            #ASDVizOptions.MomActors.spincubetmap.SetPoint2([-1.0,0.0,0.0])

            ## Calculate TCoords for texturing
            #ASDVizOptions.MomActors.spincubetmap = vtk.vtkTextureMapToSphere()
            #ASDVizOptions.MomActors.spincubetmap.SetInputConnection(ASDVizOptions.MomActors.spincube.GetOutputPort())
            ##ASDVizOptions.MomActors.spincubetmap.AutomaticSphereGenerationOff()
            ##ASDVizOptions.MomActors.spincubetmap.SetCenter([ 0.0,0.0,0.0])

            ASDVizOptions.MomActors.spincubetmap.PreventSeamOff()

            ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spincubetmap.GetOutputPort())
            ASDVizOptions.MomActors.SpinMapper.ClampingOn()
            ASDVizOptions.MomActors.SpinMapper.OrientOn()
            renWin.Render()

        if keyword=='Spheres':
            try:
                del ASDVizOptions.MomActors.spinarrow
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincube
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincones
            except:
                pass
            ASDVizOptions.MomActors.spinsphere = vtk.vtkTexturedSphereSource()
            ASDVizOptions.MomActors.spinsphere.SetRadius(0.50)
            ASDVizOptions.MomActors.spinsphere.SetThetaResolution(12)
            ASDVizOptions.MomActors.spinsphere.SetPhiResolution(12)

            # Placeholder comment for testing tangent extraction for normal textures
            #tritri = vtk.vtkTriangleFilter()
            #tritri.SetInputConnection(ASDVizOptions.MomActors.spinsphere.GetOutputPort())
            #tritan = vtk.vtkPolyDataTangents()
            #tritan.SetInputConnection(tritri.GetOutputPort())
            #ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(tritan.GetOutputPort())

            ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spinsphere.GetOutputPort())
            ASDVizOptions.MomActors.SpinMapper.ClampingOn()
            ASDVizOptions.MomActors.SpinMapper.OrientOn()
            #ASDVizOptions.MomActors.SpinMapper.OrientOff()
            renWin.Render()
        if keyword=='Arrows':
            try:
                del ASDVizOptions.MomActors.spinsphere
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincube
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincones
            except:
                pass

            # Create vectors
            ASDVizOptions.MomActors.spinarrow = vtk.vtkArrowSource()
            ASDVizOptions.MomActors.spinarrow.SetTipRadius(0.20)
            ASDVizOptions.MomActors.spinarrow.SetShaftRadius(0.10)
            ASDVizOptions.MomActors.spinarrow.SetTipResolution(12)
            ASDVizOptions.MomActors.spinarrow.SetShaftResolution(12)

            # Calculate normals for shading
            ASDVizOptions.MomActors.spinarrownormals = vtk.vtkPolyDataNormals()
            ASDVizOptions.MomActors.spinarrownormals.SetInputConnection(ASDVizOptions.MomActors.spinarrow.GetOutputPort())

            # Calculate TCoords for texturing
            ASDVizOptions.MomActors.spinarrownormalstmap = vtk.vtkTextureMapToCylinder()
            ASDVizOptions.MomActors.spinarrownormalstmap.SetInputConnection(ASDVizOptions.MomActors.spinarrownormals.GetOutputPort())
            ASDVizOptions.MomActors.spinarrownormalstmap.PreventSeamOn()

            ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spinarrownormalstmap.GetOutputPort())
            ASDVizOptions.MomActors.SpinMapper.OrientOn()
            ASDVizOptions.MomActors.SpinMapper.Update()

            renWin.Render()
        if keyword=='CenterOn':
            ASDVizOptions.MomActors.spinarrow.SetArrowOriginToCenter()
            renWin.Render()

        if keyword=='CenterOff':
            ASDVizOptions.MomActors.spinarrow.SetArrowOriginToDefault()
            renWin.Render()

        if keyword=='Cones':
            try:
                del ASDVizOptions.MomActors.spinsphere
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spincube
            except:
                pass
            try:
                del ASDVizOptions.MomActors.spinarrow
            except:
                pass

            ASDVizOptions.MomActors.spincones = vtk.vtkConeSource()
            ASDVizOptions.MomActors.spincones.SetRadius(0.50)
            ASDVizOptions.MomActors.spincones.SetHeight(1.00)
            ASDVizOptions.MomActors.spincones.SetResolution(12)

            # Calculate normals for shading
            ASDVizOptions.MomActors.spinconenormals = vtk.vtkPolyDataNormals()
            ASDVizOptions.MomActors.spinconenormals.SetInputConnection(ASDVizOptions.MomActors.spincones.GetOutputPort())

            # Calculate TCoords for texturing
            ASDVizOptions.MomActors.spinconeormalstmap = vtk.vtkTextureMapToCylinder()
            ASDVizOptions.MomActors.spinconeormalstmap.SetInputConnection(ASDVizOptions.MomActors.spinconenormals.GetOutputPort())
            ASDVizOptions.MomActors.spinconeormalstmap.PreventSeamOn()

            ASDVizOptions.MomActors.SpinMapper.SetSourceConnection(ASDVizOptions.MomActors.spinconeormalstmap.GetOutputPort())
            ASDVizOptions.MomActors.SpinMapper.OrientOn()
            ASDVizOptions.MomActors.SpinMapper.Update()
            renWin.Render()
        return
    ############################################################################
    # @brief Select the type of colormap that will be used for the different actors
    # @details Select the type of colormap that will be used for the different actors
    # It allows the user to choose between the following color schemes:
    #
    #   - Coolwarm
    #   - RdGy
    #   - Spectral
    #   - BlackBody
    # @author Jonathan Chico
    ############################################################################
    def set_colormap(self,window,flag2D,viz_type,renWin):
        """Select the type of colormap that will be used for the different actors
        It allows the user to choose between the following color schemes:
            * Coolwarm
            * RdGy
            * Spectral
            * BlackBody

        Args:
            window: QMainWindow where the visualizations are being carried out.
            flag2D: (logical) identifier indicating whether the system is in 2D or 3D.
            viz_type: (str) identifier for the different types of visualization possible in the VTK API.
            renWin: current VTK rendering window.

        Author
        ----------
        Jonathan Chico
        """
        
        import vtk

        ASDVizOptions.lut = vtk.vtkLookupTable()
        num_colors = 256
        ASDVizOptions.lut.SetNumberOfTableValues(num_colors)
        ASDVizOptions.transfer_func = vtk.vtkColorTransferFunction()
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging Coolwarm scheme by Kenneth Moreland
        #-----------------------------------------------------------------------
        if window.sender()== window.ColorMapCM and window.ColorMapCM.isChecked():
            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                ASDVizOptions.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
                ASDVizOptions.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            else:
                ASDVizOptions.transfer_func.AddRGBPoint(-1, 0.230, 0.299, 0.754)
                ASDVizOptions.transfer_func.AddRGBPoint( 1, 0.706, 0.016, 0.150)
        #-----------------------------------------------------------------------
        # Set the color to be given by the black body function
        #-----------------------------------------------------------------------
        if window.sender()== window.ColorMapBB and window.ColorMapBB.isChecked():
            ASDVizOptions.transfer_func.SetColorSpaceToRGB();
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                #ASDVizOptions.transfer_func.AddRGBPoint(0.00, 0.0, 0.0, 0.0);
                #ASDVizOptions.transfer_func.AddRGBPoint(0.25, 0.9, 0.0, 0.0);
                #ASDVizOptions.transfer_func.AddRGBPoint(0.75, 0.9, 0.9, 0.0);
                #ASDVizOptions.transfer_func.AddRGBPoint(1.00, 1.0, 1.0, 1.0);
                ASDVizOptions.transfer_func.AddRGBPoint(0.00, 0.000, 0.000, 0.000);
                ASDVizOptions.transfer_func.AddRGBPoint(0.39, 0.698, 0.133, 0.133);
                ASDVizOptions.transfer_func.AddRGBPoint(0.58, 0.890, 0.412, 0.020);
                ASDVizOptions.transfer_func.AddRGBPoint(0.89, 0.902, 0.902, 0.208);
                ASDVizOptions.transfer_func.AddRGBPoint(1.00, 1.000, 1.000, 1.000);
            else:
                ASDVizOptions.transfer_func.AddRGBPoint(-1.0, 0.0, 0.0, 0.0);
                ASDVizOptions.transfer_func.AddRGBPoint(-0.5, 0.9, 0.0, 0.0);
                ASDVizOptions.transfer_func.AddRGBPoint( 0.5, 0.9, 0.9, 0.0);
                ASDVizOptions.transfer_func.AddRGBPoint( 1.0, 1.0, 1.0, 1.0);
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging RdGy
        #-----------------------------------------------------------------------
        if window.sender()== window.ColorMapRdGy and window.ColorMapRdGy.isChecked():
            ASDVizOptions.transfer_func.SetColorSpaceToDiverging()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                ASDVizOptions.transfer_func.AddRGBPoint(0.0, 0.79216, 0.00000, 0.12549)
                ASDVizOptions.transfer_func.AddRGBPoint(0.5, 1.00000, 1.00000, 1.00000)
                ASDVizOptions.transfer_func.AddRGBPoint(1.0, 0.25098, 0.25098, 0.25098)
            else:
                ASDVizOptions.transfer_func.AddRGBPoint(-1.0, 0.79216, 0.00000, 0.12549)
                ASDVizOptions.transfer_func.AddRGBPoint( 0.0, 1.00000, 1.00000, 1.00000)
                ASDVizOptions.transfer_func.AddRGBPoint( 1.0, 0.25098, 0.25098, 0.25098)
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging spectral clor map
        #-----------------------------------------------------------------------
        if window.sender()== window.ColorMapSpectral and window.ColorMapSpectral.isChecked():
            ASDVizOptions.transfer_func.SetColorSpaceToRGB()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                ASDVizOptions.transfer_func.AddRGBPoint(0.00, 0.61961, 0.00392, 0.25882)
                ASDVizOptions.transfer_func.AddRGBPoint(0.25, 0.95686, 0.42745, 0.26275)
                ASDVizOptions.transfer_func.AddRGBPoint(0.50, 1.00000, 1.00000, 0.74902)
                ASDVizOptions.transfer_func.AddRGBPoint(0.75, 0.40000, 0.76078, 0.64706)
                ASDVizOptions.transfer_func.AddRGBPoint(1.00, 0.36863, 0.30980, 0.63529)
            else:
                ASDVizOptions.transfer_func.AddRGBPoint(-1.00, 0.61961, 0.00392, 0.25882)
                ASDVizOptions.transfer_func.AddRGBPoint(-0.50, 0.95686, 0.42745, 0.26275)
                ASDVizOptions.transfer_func.AddRGBPoint( 0.00, 1.00000, 1.00000, 0.74902)
                ASDVizOptions.transfer_func.AddRGBPoint( 0.50, 0.40000, 0.76078, 0.64706)
                ASDVizOptions.transfer_func.AddRGBPoint( 1.00, 0.36863, 0.30980, 0.63529)
        #-----------------------------------------------------------------------
        # Construct the lut with the selected colomap
        #-----------------------------------------------------------------------
        for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
            cc = ASDVizOptions.transfer_func.GetColor(ss)
            ASDVizOptions.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
        ASDVizOptions.lut.SetTableRange(0.0,1.0)
        ASDVizOptions.lut.Build()
        #-----------------------------------------------------------------------
        # Color the actors depending of the type of visualization
        #-----------------------------------------------------------------------
        if viz_type=='M':
            if (flag2D and viz_type=='M') or (flag2D and viz_type=='E') or viz_type=='N':
                ASDVizOptions.MomActors.MagDensMap.SetLookupTable(ASDVizOptions.lut)
                ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.lut)
                ASDVizOptions.GenActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
                ASDVizOptions.GenActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
            else:
                ASDVizOptions.MomActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
                ASDVizOptions.MomActors.SpinMapper.SetLookupTable(ASDVizOptions.transfer_func)
                ASDVizOptions.GenActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
                ASDVizOptions.GenActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)
        elif viz_type=='N':
            ASDVizOptions.NeighActors.NeighMapper.SetLookupTable(ASDVizOptions.lut)
            ASDVizOptions.GenActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
            ASDVizOptions.GenActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
        elif viz_type=='E':
            if flag2D:
                ASDVizOptions.EneActors.EneDensMap.SetLookupTable(ASDVizOptions.lut)
                ASDVizOptions.GenActors.scalar_bar.SetLookupTable(ASDVizOptions.lut)
                ASDVizOptions.GenActors.clipperMapper.SetLookupTable(ASDVizOptions.lut)
            else:
                ASDVizOptions.EneActors.volumeProperty.SetColor(ASDVizOptions.transfer_func)
                ASDVizOptions.GenActors.scalar_bar.SetLookupTable(ASDVizOptions.transfer_func)
                ASDVizOptions.GenActors.clipperMapper.SetLookupTable(ASDVizOptions.transfer_func)
        #-----------------------------------------------------------------------
        # Render the scene
        #-----------------------------------------------------------------------
        renWin.Render()
        return

    ############################################################################
    # @brief Select the type of colormap that will be used for the different actors
    # @details Select the type of colormap that will be used for the different actors
    # It allows the user to choose between the following color schemes:
    #
    #   - Coolwarm
    #   - RdGy
    #   - Spectral
    #   - BlackBody
    # @author Jonathan Chico
    ############################################################################
    def set_colormap_db(self,mapnum,window,flag2D,viz_type,renWin):
        """Select the type of colormap that will be used for the different actors
        It allows the user to choose between the following color schemes:
            * Coolwarm (mapnum 0)
            * RdGy (mapnum 1)
            * Spectral (mapnum 2)
            * BlackBody (mapnum 3)

        Args:
            window: QMainWindow where the visualizations are being carried out.
            flag2D: (logical) identifier indicating whether the system is in 2D or 3D.
            viz_type: (str) identifier for the different types of visualization possible in the VTK API.
            renWin: current VTK rendering window.

        Author
        ----------
        Jonathan Chico
        """
        
        import vtk

        #self.lut = vtk.vtkLookupTable()
        num_colors = 256
        self.lut.SetNumberOfTableValues(num_colors)
        self.transfer_func = vtk.vtkColorTransferFunction()
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging Coolwarm scheme by Kenneth Moreland
        #-----------------------------------------------------------------------
        if mapnum == 0:
            self.transfer_func.SetColorSpaceToDiverging()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                self.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
                self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            else:
                self.transfer_func.AddRGBPoint(-1, 0.230, 0.299, 0.754)
                self.transfer_func.AddRGBPoint( 1, 0.706, 0.016, 0.150)
        #-----------------------------------------------------------------------
        # Set the color to be given by the black body function
        #-----------------------------------------------------------------------
        if mapnum == 1:
            self.transfer_func.SetColorSpaceToRGB();
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                self.transfer_func.AddRGBPoint(0.00, 0.000, 0.000, 0.000);
                self.transfer_func.AddRGBPoint(0.39, 0.698, 0.133, 0.133);
                self.transfer_func.AddRGBPoint(0.58, 0.890, 0.412, 0.020);
                self.transfer_func.AddRGBPoint(0.89, 0.902, 0.902, 0.208);
                self.transfer_func.AddRGBPoint(1.00, 1.000, 1.000, 1.000);
            else:
                self.transfer_func.AddRGBPoint(-1.0, 0.0, 0.0, 0.0);
                self.transfer_func.AddRGBPoint(-0.5, 0.9, 0.0, 0.0);
                self.transfer_func.AddRGBPoint( 0.5, 0.9, 0.9, 0.0);
                self.transfer_func.AddRGBPoint( 1.0, 1.0, 1.0, 1.0);
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging RdGy
        #-----------------------------------------------------------------------
        if mapnum == 2:
            self.transfer_func.SetColorSpaceToDiverging()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                self.transfer_func.AddRGBPoint(0.0, 0.79216, 0.00000, 0.12549)
                self.transfer_func.AddRGBPoint(0.5, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 0.25098, 0.25098, 0.25098)
            else:
                self.transfer_func.AddRGBPoint(-1.0, 0.79216, 0.00000, 0.12549)
                self.transfer_func.AddRGBPoint( 0.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint( 1.0, 0.25098, 0.25098, 0.25098)
        #-----------------------------------------------------------------------
        # Set the color map to be given by the diverging spectral clor map
        #-----------------------------------------------------------------------
        if mapnum == 3:
            self.transfer_func.SetColorSpaceToRGB()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                self.transfer_func.AddRGBPoint(0.00, 0.61961, 0.00392, 0.25882)
                self.transfer_func.AddRGBPoint(0.25, 0.95686, 0.42745, 0.26275)
                self.transfer_func.AddRGBPoint(0.50, 1.00000, 1.00000, 0.74902)
                self.transfer_func.AddRGBPoint(0.75, 0.40000, 0.76078, 0.64706)
                self.transfer_func.AddRGBPoint(1.00, 0.36863, 0.30980, 0.63529)
            else:
                self.transfer_func.AddRGBPoint(-1.00, 0.61961, 0.00392, 0.25882)
                self.transfer_func.AddRGBPoint(-0.50, 0.95686, 0.42745, 0.26275)
                self.transfer_func.AddRGBPoint( 0.00, 1.00000, 1.00000, 0.74902)
                self.transfer_func.AddRGBPoint( 0.50, 0.40000, 0.76078, 0.64706)
                self.transfer_func.AddRGBPoint( 1.00, 0.36863, 0.30980, 0.63529)
        #-----------------------------------------------------------------------
        # High-jacking this scheme for Single colors
        if mapnum == -1:
            self.transfer_func.SetColorSpaceToRGB()
            if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
                self.transfer_func.AddRGBPoint(0.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 1.00000, 1.00000, 1.00000)
            else:
                self.transfer_func.AddRGBPoint(-1.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint( 1.0, 1.00000, 1.00000, 1.00000)
        # Construct the lut with the selected colomap
        #-----------------------------------------------------------------------
        for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
            cc = self.transfer_func.GetColor(ss)
            self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
        self.lut.Build()
        return

    ############################################################################
    # @brief Set the background color given RGB values
    # @author Anders Bergman
    ############################################################################
    def set_RGBbackground(self,rgb,ren,renWin):

        nrgb = [ i/255.0 for i in rgb ]
        ren.SetBackground(nrgb)


        return
    ############################################################################
    # @brief Setup a single color colormap 
    # @details Used since glyphs are colored after colormaps and not single colors
    # @author Anders Bergman
    ############################################################################
    def set_RGBcolor(self,rgb,window,flag2D,viz_type,renWin):
        
        import vtk


        num_colors = 2
        red = rgb[0]/255.0
        green = rgb[1]/255.0
        blue = rgb[2]/255.0

        self.lut.SetNumberOfTableValues(num_colors)
        self.transfer_func = vtk.vtkColorTransferFunction()
        self.transfer_func.SetColorSpaceToRGB();
        if (viz_type=='M') or (viz_type=='E') or viz_type=='N':
            self.transfer_func.AddRGBPoint(0.00, red, green, blue)
            self.transfer_func.AddRGBPoint(1.00, red, green, blue)
        else:
            self.transfer_func.AddRGBPoint(-1.0, red, green, blue)
            self.transfer_func.AddRGBPoint( 1.0, red, green, blue)

        for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
            cc = self.transfer_func.GetColor(ss)
            self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)

        self.lut.Build()
        return