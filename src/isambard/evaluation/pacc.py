"""This module contains tools for the Parametric Analysis of Coiled Coils."""

import numpy

from ampal.analyse_protein import reference_axis_from_chains, alpha_angles, crick_angles,\
    polymer_to_reference_axis_distances, polypeptide_vector
from ampal.pseudo_atoms import Primitive
from ampal.geometry import is_acute

class PACCAnalysis(object):
    def __init__(self, coiledcoil):
        """Class for the parametric analysis of coiled coils.

        Currently only functions for parallel and ap blunt ended assemblies.

        Parameters
        ----------
        coiledcoil: Assembly
            Must contain only the coiled coil polypeptides which need to be of
            equal length.
        """
        len_set = set([len(x) for x in coiledcoil])
        if len(len_set) != 1:
            raise ValueError('The helices of the coiled coil must be of equal length.')
        self.cc_len = len_set.pop()
        self.cc = coiledcoil
        self.ra = reference_axis_from_chains(self.cc)
        # create flipped axis
        self.ra_flipped = Primitive.from_coordinates(numpy.flipud(self.ra.coordinates))
        ref_polypeptide_vec = polypeptide_vector(self.cc[0])

        for ch in self.cc:
            ch_polypeptide_vec = polypeptide_vector(ch)
            # if both vectors point in the same direction (angle less than 90 deg)
            if is_acute(ref_polypeptide_vec, ch_polypeptide_vec):
                ref_ax = self.ra
            else:
                ref_ax = self.ra_flipped
            polymer_to_reference_axis_distances(ch, ref_ax)
            alpha_angles(ch, ref_ax)
            crick_angles(ch, ref_ax)

        self.radii_layers = []
        self.alpha_layers = []
        self.ca_layers = []
        self.gather_layer_info()

    def gather_layer_info(self):
        """Extracts the tagged coiled-coil parameters for each layer."""
        for i in range(len(self.cc[0])):
            layer_radii = [x[i].tags['distance_to_ref_axis'] for x in self.cc]
            self.radii_layers.append(layer_radii)
            layer_alpha = [x[i].tags['alpha_angle_ref_axis'] for x in self.cc]
            self.alpha_layers.append(layer_alpha)
            layer_ca = [x[i].tags['crick_angle_ref_axis'] for x in self.cc]
            self.ca_layers.append(layer_ca)
        return

    @staticmethod
    def calc_average_parameters(parameter_layers):
        """Takes a group of equal length lists and averages them across each index.

        Returns
        -------
        mean_layers: [float]
            List of values averaged by index
        overall_mean: float
            Mean of the averaged values.
        """
        mean_layers = [numpy.mean(x) if x[0] else 0 for x in parameter_layers]
        overall_mean = numpy.mean([x for x in mean_layers if x])
        return mean_layers, overall_mean

    def heptad_register(self):
        """Returns the calculated register of the coiled coil and the fit quality."""
        base_reg = 'abcdefg'
        exp_base = base_reg * (self.cc_len//7+2)
        ave_ca_layers = self.calc_average_parameters(self.ca_layers)[0][:-1]
        reg_fit = fit_heptad_register(ave_ca_layers)
        hep_pos = reg_fit[0][0]
        return exp_base[hep_pos:hep_pos+self.cc_len], reg_fit[0][1:]

    def generate_report(self):
        """Generates a report on the coiled coil parameters.

        Returns
        -------
        report: str
            A string detailing the register and parameters of the coiled coil.
        """
        # Find register
        lines = ['Register Assignment\n-------------------']
        register, fit = self.heptad_register()
        lines.append('{}\n{}\n'.format(register, '\n'.join(self.cc.sequences)))
        lines.append('Fit Quality - Mean Angular Discrepancy = {:3.2f} (Std Dev = {:3.2f})\n'.format(*fit))
        # Find coiled coil parameters
        lines.append('Coiled Coil Parameters\n----------------------')
        layer_info = (self.radii_layers, self.alpha_layers, self.ca_layers)
        r_layer_aves, a_layer_aves, c_layer_aves = [self.calc_average_parameters(x) for x in layer_info]
        start_line = ['Res#'.rjust(5), 'Radius'.rjust(9), 'Alpha'.rjust(9), 'CrAngle'.rjust(9)]
        lines.append(''.join(start_line))
        for i in range(len(r_layer_aves[0])):
            residue = '{:>5}'.format(i+1)
            average_r = '{:+3.3f}'.format(r_layer_aves[0][i]).rjust(9)
            average_a = '{:+3.3f}'.format(a_layer_aves[0][i]).rjust(9)
            average_c = '{:+3.3f}'.format(c_layer_aves[0][i]).rjust(9)
            line = [residue, average_r, average_a, average_c]
            lines.append(''.join(line))
        # Average for assembly
        lines.append('-'*32)
        residue = '  Ave'
        average_r = '{:+3.3f}'.format(r_layer_aves[1]).rjust(9)
        average_a = '{:+3.3f}'.format(a_layer_aves[1]).rjust(9)
        average_c = '{:+3.3f}'.format(c_layer_aves[1]).rjust(9)
        line = [residue, average_r, average_a, average_c]
        lines.append(''.join(line))
        # Std dev
        residue = 'Std D'
        std_d_r = '{:+3.3f}'.format(numpy.std(r_layer_aves[0])).rjust(9)
        std_d_a = '{:+3.3f}'.format(numpy.std(a_layer_aves[0][:-1])).rjust(9)
        std_d_c = '{:+3.3f}'.format(numpy.std(c_layer_aves[0][:-1])).rjust(9)
        line = [residue, std_d_r, std_d_a, std_d_c]
        lines.append(''.join(line))
        return '\n'.join(lines)


def fit_heptad_register(crangles):
    """Attempts to fit a heptad repeat to a set of Crick angles.

    Parameters
    ----------
    crangles: [float]
        A list of average Crick angles for the coiled coil.

    Returns
    -------
    fit_data: [(float, float, float)]
        Sorted list of fits for each heptad position.
    """
    crangles = [x if x > 0 else 360 + x for x in crangles]
    hept_p = [x * (360.0 / 7.0) + ((360.0 / 7.0) / 2.0) for x in range(7)]
    ideal_crangs = [
        hept_p[0],
        hept_p[2],
        hept_p[4],
        hept_p[6],
        hept_p[1],
        hept_p[3],
        hept_p[5]
    ]
    full_hept = len(crangles) // 7
    ideal_crang_list = ideal_crangs * (full_hept + 2)  # This is dirty, too long but trimmed with zip
    fitting = []
    for i in range(7):
        ang_pairs = zip(crangles, ideal_crang_list[i:])
        ang_diffs = [abs(y - x) for x, y in ang_pairs]
        fitting.append((i, numpy.mean(ang_diffs), numpy.std(ang_diffs)))
    return sorted(fitting, key=lambda x: x[1])


__author__ = 'Christopher W. Wood'
