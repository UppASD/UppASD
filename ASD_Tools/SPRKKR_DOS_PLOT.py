class SPRKKRPlotDos():
    """Script to plot the DOS and orbital projected DOS obtained from SPRKKR.
    It works by reading the .dos raw data file created by SPRKKR versions
    6.3, 7.1 and 7.7 then it transforms the data to appropriate arrays which
    then are plotted.
    The script generates legends with appropriate names, allows one to choose
    the energy range of plotting, transforms the energy scale to eVs and allows
    one to decide whether the DOS of an element should be plotted.
    As an input it will only ask for the name of the .dos file that must be
    plotted.
    Tested with python 2.7 and 3.6
    """

    def __init__(self):
        self.fig_size = (16, 10)
        self.fig_dpi = 800
        self.fig_font = 28
        self.ax_font = 20
        self.ene_range = [-7.5, 7.5]
        self.exclude_list = ['Vc']
        self.SPRKKR_char_len = 10
        self.DOS_units = 'ev'
        self.DOS_data = []
        self.fermi_energy = []
        self.atom_types = []
        self.element_bare = []
        self.num_types = 1
        self.sys_name = None
        self.fig_kwargs = None
        self.fig_type = 'filled'
        return

    def SPRKKR_DOS_setup(self):
        self.get_args()
        self.set_default_fig_dict()
        self.get_extras()
        return

    def SPRKKR_DOS_driver(self):
        self.SPRKKR_DOS_parse()
        self.SPRKKR_DOS_plot()
        return

    def ry_to_ev(self):
        from scipy import constants
        return constants.value(u'Rydberg constant times hc in eV')

    def is_latex(self):
        from distutils.spawn import find_executable
        return find_executable('latex')

    def set_default_fig_dict(self):
        if self.fig_type.lower() == 'filled':
            self.fig_kwargs = {'alpha': 0.5, 'linewidth': 3, 'linestyle': '-'}
        else:
            self.fig_kwargs = {'alpha': 1.0, 'linewidth': 3, 'linestyle': '-'}
        return

    # Functions to get the arguments from the CLI and pass them to the program
    def get_args(self):
        """Function to handle the definition of runtime  variables which can be
        passed on call.
        """
        import argparse
        import glob
        # Setting the information displayed when the --help or -h is displayed
        _msg = 'Plots SPRKKR DOS files.'
        parser = argparse.ArgumentParser(description=_msg)
        # input file argument
        _msg = 'Name of the file containing the input parameters'
        parser.add_argument("--extras", "-x", help=_msg, default='DOS_inp.yml',
                            metavar="FILE", type=str)
        # output file argument
        _msg = 'Name of the SPRKKR .dos file'
        parser.add_argument("--input", "-i", help=_msg,
                            default=glob.glob('*DOS.dos')[0], metavar="FILE",
                            type=str)
        # log file argument
        _msg = 'Name of the log file'
        parser.add_argument("--log", "-l", help=_msg, default='cell_file.log',
                            metavar="FILE", type=str)
        # error file argument
        _msg = 'Name of the error file'
        parser.add_argument("--error", "-e", help=_msg,
                            default='cell_file.err', metavar="FILE", type=str)
        # verbose argument
        _msg = 'Print messages to screen.'
        parser.add_argument('--verbose', "-v", help=_msg, default=False,
                            action="store_true")
        # Label for the full plot
        _msg = 'Adds a label for the full plot.'
        parser.add_argument('--flabel', "-fl", help=_msg, default=None,
                            type=str)
        # Label for the full plot
        _msg = 'Adds a label for the types plot.'
        parser.add_argument('--tlabel', "-tl", help=_msg, default=None,
                            type=str)
        # Parse the arguments so that they can be used
        self.args = parser.parse_args()
        return

    def get_energies(self):
        if self.DOS_units.lower() == 'ev':
            return self.DOS_data[0]*self.ry_to_ev()
        else:
            return self.DOS_data[0]

    def get_full_DOS(self):
        if self.DOS_units.lower() == 'ev':
            return self.DOS_data[1]*(1.0/self.ry_to_ev())
        else:
            return self.DOS_data[1]

    def get_atom_DOS(self):
        if self.DOS_units.lower() == 'ev':
            return self.DOS_data[2]*(1.0/self.ry_to_ev())
        else:
            return self.DOS_data[2]

    def get_fermi_energy(self):
        if self.DOS_units.lower() == 'ev':
            return self.fermi_energy*self.ry_to_ev()
        else:
            return self.fermi_energy

    def get_atom_types(self):
        return self.atom_types

    # Create a function that reads the lattice vector and basis
    def get_extras(self):
        """ Function to read the input parameters needed for the plotting of
        the DOS from SPRKKR.
        """
        import sys
        import yaml

        try:
            with open(self.args.extras, 'r') as stream:
                dos_extras = yaml.load(stream, Loader=yaml.FullLoader)
            # Function to transform all the dictionary functions to the
            # appropriate types needed for the rest of the functions
            self.init_var(dos_extras)

        except FileNotFoundError:
            pass
        except IOError:
            sys.exit("Problem reading " + self.args.extras +
                     " check the format")
        return

    def init_var(self, dos_extras):
        import numpy as np
        if 'Fig' in dos_extras:
            if 'dpi' in dos_extras['Fig']:
                self.fig_dpi = int(dos_extras['Fig']['dpi'])
            if 'figsize' in dos_extras['Fig']:
                self.fig_size = dos_extras['Fig']['figsize']
            if 'fontsize' in dos_extras['Fig']:
                self.fig_font = float(dos_extras['Fig']['fontsize'])
            if 'linewidth' in dos_extras['Fig']:
                self.fig_kwargs['linewidth'] = \
                        int(dos_extras['Fig']['linewidth'])
            if 'alpha' in dos_extras['Fig']:
                self.fig_kwargs['alpha'] = float(dos_extras['Fig']['alpha'])
            if 'linestyle' in dos_extras['Fig']:
                self.fig_kwargs['linestyle'] = \
                        str(dos_extras['Fig']['linestyle'])
            if 'type' in dos_extras['Fig']:
                self.fig_type = str(dos_extras['Fig']['type']).lower()
            if 'full_label' in dos_extras['Fig']\
               and self.args.flabel is not None:
                self.args.flabel = dos_extras['Fig']['full_label']
        if 'Data' in dos_extras:
            if 'ene_bounds' in dos_extras['Data']:
                self.ene_range = np.asarray(dos_extras['Data']['ene_bounds'],
                                            dtype=np.float64)
            if 'exclude' in dos_extras['Data']:
                self.exclude_list = dos_extras['Data']['exclude']
            if 'units' in dos_extras['Data']:
                self.DOS_units = str(dos_extras['Data']['units'])
        if 'Misc' in dos_extras:
            if 'char_len' in dos_extras['Misc']:
                self.SPRKKR_char_len = int(dos_extras['Misc']['char_len'])
        return

    def SPRKKR_DOS_parse(self):
        """Function wrapper for reading and parsing the DOS file produced by
        SPRKKR
        """
        # This is the name of the file that we want to plot
        self.DOS_file = open(self.args.input)
        # Paring the header which contains relevant information of the system
        self.SPRKKR_DOS_parse_header()
        # Parsing the actual data
        self.SPRKKR_DOS_parse_data()
        return

    def SPRKKR_DOS_parse_header(self):
        #######################################################################
        # Initialization of certain variables
        #######################################################################
        # Set a flag to check the start of a new type of data to be False
        _data_found = False
        # Initiate the counter
        count = 0
        #######################################################################
        # Check for certain key data in the file
        #######################################################################
        while not _data_found:
            count = count + 1
            _line = self.DOS_file.readline()
            _data = str.split(_line)
            # Find the name of the system to be used in the output name
            if len(_data) > 0 and _data[0] == 'SYSTEM':
                self.sys_name = str(_data[1])
            # Find the number of non-equivalent atoms present in the system
            if len(_data) > 0 and _data[0] == 'NT_eff':
                self.num_types = int(_data[1])
            # Find the number of sites present in the system
            if len(_data) > 0 and _data[0] == 'NQ_eff':
                self.num_chem = int(_data[1])
            # Find the number of energy points found in the calculation
            if len(_data) > 0 and _data[0] == 'NE':
                self.num_ene = int(_data[1])
            # Find the Fermi energy of the system
            if len(_data) > 0 and _data[0] == 'EFERMI':
                self.fermi_energy = float(_data[1])
            # Find the angular momentum expansion
            if len(_data) > 0 and _data[0] == 'IQ':
                for ii in range(0, self.num_chem):
                    _line = self.DOS_file.readline()
                    _data = str.split(_line)
                    self.lmax = int(_data[1])
            # Find the labels associated with the different atoms types present
            # in the system
            if len(_data) > 0 and _data[0] == 'IT':
                ii = 0
                while ii < self.num_types:
                    _line = self.DOS_file.readline()
                    _data = str.split(_line)
                    if len(_data) > 2:
                        curr_string = "$" + str(_data[1]) + "$"
                        ind = str(_data[1]).find('_')
                        if ind > 0:
                            self.element_bare.append(str(_data[1][0:ind]))
                        else:
                            self.element_bare.append(str(_data[1][0:]))
                        self.atom_types.append(curr_string)
                        ii += 1
                _data_found = True
        #######################################################################
        # Rewind the file
        #######################################################################
        self.DOS_file.seek(0)
        return

    def SPRKKR_DOS_parse_data(self):
        import numpy as np
        import itertools
        import matplotlib.pyplot as plt
        #######################################################################
        # Set the discriminating flag to False again
        #######################################################################
        _data_found = False
        #######################################################################
        # Check when the actual data starts
        #######################################################################
        while not _data_found:
            line = self.DOS_file.readline()
            _data = str.split(line)
            if len(_data) > 0 and _data[0] == 'DOS-FMT:':
                _data_found = True
        #######################################################################
        # The reading of the actual DOS data begins
        #######################################################################
        _data = []  # Array that contains the actual full data
        #######################################################################
        # Read line by line splitting the data in accordance with what is
        # needed
        #######################################################################
        _ene = []
        _data_full = []
        count = 0
        while line:
            line = self.DOS_file.readline()
            if len(line) > 0:
                if line[0] != ' ':
                    count = 1
                    line = line.replace(" ", "")
                    curr_data =\
                        [float(line[ii: ii + self.SPRKKR_char_len])
                               for ii in range(0, len(line)-1,
                                               self.SPRKKR_char_len)]
                else:
                    tmp = self.SPRKKR_char_len*str(0)\
                            + line[self.SPRKKR_char_len:]
                    tmp = tmp.replace(" ", "")
                    curr_data = [float(tmp[ii: ii + self.SPRKKR_char_len])
                                 for ii in range(0, len(tmp)-1,
                                                 self.SPRKKR_char_len)]
                    count = count + 1
                # Append the data to an actually usuable shape
                _data.append(curr_data)
        self.DOS_file.close()
        #######################################################################
        # This ensures that the full range of data associated with a given
        # energy value is in a single array for convenience
        #######################################################################
        for ii in range(0, len(_data)-1):
            if ii % count == 0:
                _curr = _data[ii][1: len(_data[0])]
                for kk in range(1, count):
                    for jj in range(1, len(_data[ii + kk][:])):
                        _curr = np.hstack((_curr, _data[ii + kk][jj]))
                _data_full.append(_curr)
                _ene.append((float(_data[ii][0]) - self.fermi_energy))
        self.DOS_data.append(np.asarray(_ene))
        self.DOS_data.append(np.asarray(_data_full))
        #######################################################################
        # Creates an atom resolved array with the following dimensions
        # (Number of energies, Number of atoms types, (s,p,d,f,Total), Spin)
        #######################################################################
        _atm_res = np.zeros([self.num_ene, self.num_types, self.lmax + 1, 2])
        # Stores the data in the site resolved array for easy plotting
        for ene, nt, spin, lmn in itertools.product(range(0, self.num_ene),
                                                    range(0, self.num_types),
                                                    range(0, 2),
                                                    range(0, self.lmax)):
            column_num = int(lmn + self.lmax*spin + self.lmax*2*nt + 1)
            _atm_res[ene, nt, lmn, spin] = float(_data_full[ene][column_num])
        _atm_res[:, :, self.lmax, :] = np.sum(_atm_res[:, :, :, :], axis=2)
        self.DOS_data.append(_atm_res)
        del _ene, _data_full, _curr, _atm_res, _data
        return

    def get_num_colors(self):
        """Finds the number of colors needed for the plot after one takes care
        of the exclusion list
        """
        return len([x for x in self.element_bare
                    if x not in self.exclude_list])

    def get_DOS_max_type(self, ene_indx=-1):
        import numpy as np
        """Finds the maximum DOS per atom type
        """
        return [np.maximum(self.get_atom_DOS()
                           [ene_indx, nt, self.lmax, 0].max(),
                           self.get_atom_DOS()
                           [ene_indx, nt, self.lmax, 1].max())
                for nt in range(0, self.num_types)]

    def SPRKKR_DOS_plot(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import cm as cm

        if self.is_latex():
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

        if self.DOS_units.lower() == 'ev':
            x_label = r'E-E$_F$ [eV]'
            y_label = [r'n$_\downarrow$ [sts./eV]', r'n$_\uparrow$ [sts./eV]']
        else:
            x_label = r'E-E$_F$ [Ry]'
            y_label = [r'n$_\downarrow$ [sts./Ry]', r'n$_\uparrow$ [sts./Ry]']

        # Find the energy range
        ene_indx = np.where((self.get_energies() >= self.ene_range[0]) &
                            (self.get_energies() < self.ene_range[1]))[0]

        num_colors = self.get_num_colors()
        max_DOS_NT = self.get_DOS_max_type(ene_indx)
        max_DOS = np.max(max_DOS_NT)
        colors = cm.Paired(np.linspace(0, 1, num_colors))
        fig, (ax1, ax2) = \
            plt.subplots(2, 1, sharex=True, figsize=self.fig_size)
        #######################################################################
        # If one chooses to one can eliminate the Vc DOS from the total plot
        #######################################################################
        prune_atom_types = []
        count = 0
        for nt in range(0, self.num_types):
            if self.element_bare[nt] not in self.exclude_list:
                ###############################################################
                # Plots the top panel
                ###############################################################
                if self.fig_type.lower() == 'filled':
                    ax1.fill_between(self.get_energies()[ene_indx],
                                     self.get_atom_DOS()
                                     [ene_indx, nt, self.lmax, 0],
                                     facecolor=colors[count],
                                     edgecolor=colors[count],
                                     **self.fig_kwargs)
                    ###########################################################
                    # Plots the lower panel
                    ###########################################################
                    ax2.fill_between(self.get_energies()[ene_indx],
                                     self.get_atom_DOS()
                                     [ene_indx, nt, self.lmax, 1],
                                     facecolor=colors[count],
                                     edgecolor=colors[count],
                                     **self.fig_kwargs)
                else:
                    ax1.plot(self.get_energies()[ene_indx],
                             self.get_atom_DOS()[ene_indx, nt, self.lmax, 0],
                             color=colors[count], **self.fig_kwargs)
                    ###########################################################
                    # Plots the lower panel
                    ###########################################################
                    ax2.plot(self.get_energies()[ene_indx],
                             self.get_atom_DOS()[ene_indx, nt, self.lmax, 1],
                             color=colors[count], **self.fig_kwargs)
                prune_atom_types.append(self.get_atom_types()[nt])
                count += 1
        #######################################################################
        # Creates the legend that will be displayed
        #######################################################################
        ax1.legend(prune_atom_types, fontsize=self.ax_font, loc='upper right')
        #######################################################################
        # Plotting options
        #######################################################################
        # Set the limits in the y-axis of the plot
        ax1.set_ylim([0.00, max_DOS*1.1])
        # Sets the limits in the x-axis of the plot
        ax1.set_xlim([self.get_energies()[ene_indx].min(),
                      self.get_energies()[ene_indx].max()])
        ax1.set_ylabel(y_label[0], fontsize=self.fig_font)
        ax1.axvline(0, color='black', linestyle='--', linewidth=1.25)
        ax1.tick_params(axis='x', colors='black', labelsize=self.ax_font,
                        width=0)
        ax1.tick_params(axis='y', colors='black', labelsize=self.ax_font,
                        width=2)
        for axis in ['top', 'left', 'right']:
            ax1.spines[axis].set_linewidth(3)
        for axis in ['bottom']:
            ax1.spines[axis].set_linewidth(0)
        # Set the limits in the y-axis of the plot
        ax2.set_ylim([0.00, max_DOS*1.1])
        # Sets the limits in the x-axis of the plot
        ax2.set_xlim([self.get_energies()[ene_indx].min(),
                      self.get_energies()[ene_indx].max()])
        # Inverts the y-axis to obtain a traditional DOS style plot
        ax2.invert_yaxis()
        ax2.axvline(0, color='black', linestyle='--', linewidth=1.25)
        ax2.tick_params(axis='x', colors='black', labelsize=self.ax_font,
                        width=2)
        ax2.tick_params(axis='y', colors='black', labelsize=self.ax_font,
                        width=2)
        ax2.set_xlabel(x_label, fontsize=self.fig_font)
        ax2.set_ylabel(y_label[1], fontsize=self.fig_font)
        fig.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        for axis in ['bottom', 'left', 'right']:
            ax2.spines[axis].set_linewidth(3)
        for axis in ['top']:
            ax2.spines[axis].set_linewidth(0)
        # Add a label to the figure
        if self.args.flabel is not None:
            if len(self.args.flabel) > 0:
                plt.text(0.075, 0.85, self.args.flabel,
                         fontsize=self.fig_font,
                         transform=plt.gcf().transFigure, weight='bold')
        # Set the size of the plot
        # Define the current filename
        fig_name = 'Total_DOS_' + self.sys_name + '.pdf'
        # Save the figure to file
        fig.savefig(fig_name, transparent=False, dpi=self.fig_dpi,
                    bbox_inches='tight')
        plt.clf()
        #######################################################################
        # Create an orbital resolved DOS plot for each atom type
        #######################################################################
        # Define an array with the name of the orbitals
        orbitals = ['s', 'p', 'd', 'f']
        # Do a loop over all the atom types present in the system
        colors = cm.Paired(np.linspace(0, 1, self.lmax))
        for nt in range(0, self.num_types):
            if self.element_bare[nt] not in self.exclude_list:
                loc_max_site = np.zeros(self.lmax, dtype=np.float64)
                # Find the maximum value of the DOS for current atom type
                for lnm in range(0, self.lmax):
                    loc_max_site[lnm] = \
                        np.maximum(self.get_atom_DOS()
                                   [ene_indx, nt, lnm, 0].max(),
                                   self.get_atom_DOS()
                                   [ene_indx, nt, lnm, 1].max())
                loc_max = np.max(loc_max_site)
                fig, (ax1, ax2) =\
                    plt.subplots(2, 1, sharex=True, figsize=self.fig_size)
                # Plots the top panel of the total DOS for all the orbitals
                if self.fig_type.lower() == 'filled':
                    for ii in range(0, self.lmax):
                        ax1.fill_between(self.get_energies()[ene_indx],
                                         self.get_atom_DOS()
                                         [ene_indx, nt, ii, 0],
                                         facecolor=colors[ii],
                                         edgecolor=colors[ii],
                                         **self.fig_kwargs)
                else:
                    for ii in range(0, self.lmax):
                        ax1.plot(self.get_energies()[ene_indx],
                                 self.get_atom_DOS()[ene_indx, nt, ii, 0],
                                 color=colors[ii], **self.fig_kwargs)
                ax1.legend(orbitals, fontsize=self.fig_font, loc='upper right')
                # Sets the range of the y-axis
                ax1.set_ylim([0, loc_max*1.1])
                # Sets the limits of the x-axi
                ax1.set_xlim([self.get_energies()[ene_indx].min(),
                              self.get_energies()[ene_indx].max()])
                ax1.set_ylabel(y_label[0], fontsize=self.fig_font)
                ax1.axvline(0, color='black', linestyle='--', linewidth=1.25)
                ax1.tick_params(axis='x', colors='black',
                                labelsize=self.ax_font, width=0)
                ax1.tick_params(axis='y', colors='black',
                                labelsize=self.ax_font, width=2)
                for axis in ['top', 'left', 'right']:
                    ax1.spines[axis].set_linewidth(3)
                for axis in ['bottom']:
                    ax1.spines[axis].set_linewidth(0)
                # Plots the bottom panel of the total DOS for all the orbitals
                if self.fig_type.lower() == 'filled':
                    for ii in range(0, self.lmax):
                        ax2.fill_between(self.get_energies()[ene_indx],
                                         self.get_atom_DOS()
                                         [ene_indx, nt, ii, 1],
                                         facecolor=colors[ii],
                                         edgecolor=colors[ii],
                                         **self.fig_kwargs)
                else:
                    for ii in range(0, self.lmax):
                        ax2.plot(self.get_energies()[ene_indx],
                                 self.get_atom_DOS()[ene_indx, nt, ii, 1],
                                 color=colors[ii], **self.fig_kwargs)
                # Set the limits in the y-axis of the plot
                ax2.set_ylim([0, loc_max*1.1])
                # Sets the limits in the x-axis of the plot
                ax2.set_xlim([self.get_energies()[ene_indx].min(),
                              self.get_energies()[ene_indx].max()])
                # Inverts the y-axis to obtain a traditional DOS style plot
                ax2.invert_yaxis()
                ax2.axvline(0, color='black', linestyle='--', linewidth=1.25)
                ax2.tick_params(axis='x', colors='black',
                                labelsize=self.ax_font, width=2)
                ax2.tick_params(axis='y', colors='black',
                                labelsize=self.ax_font, width=2)
                ax2.set_xlabel(x_label, fontsize=self.fig_font)
                ax2.set_ylabel(y_label[1], fontsize=self.fig_font)
                fig.subplots_adjust(hspace=0)
                plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],
                         visible=False)
                for axis in ['bottom', 'left', 'right']:
                    ax2.spines[axis].set_linewidth(3)
                for axis in ['top']:
                    ax2.spines[axis].set_linewidth(0)
                # Define the current filename
                _nt = len(self.get_atom_types()[nt])-1
                # Add the labels to the plots
                if self.args.tlabel is not None:
                    if len(self.args.tlabel) > 0:
                        plt.text(0.075, 0.85, self.args.tlabel,
                                 fontsize=self.fig_font,
                                 transform=plt.gcf().transFigure,
                                 weight='bold')
                fig_name = 'Local_DOS_' + self.sys_name + '_' +\
                    self.get_atom_types()[nt][1: _nt] + '.pdf'
                # Save the figure to file
                fig.savefig(fig_name, transparent=False, dpi=self.fig_dpi,
                            bbox_inches='tight')
                plt.clf()
        return


def main():

    SPRKKRDOS = SPRKKRPlotDos()
    SPRKKRDOS.SPRKKR_DOS_setup()
    SPRKKRDOS.SPRKKR_DOS_driver()

    return


if __name__ == '__main__':
    main()
