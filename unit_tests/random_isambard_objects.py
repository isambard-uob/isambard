import random
import numpy

import isambard_dev as isambard


def random_angles(n=1, min_val=0, max_val=180, radians=False):
    angles = [(random.random() * random.choice(range(abs(max_val - min_val)))) + min_val for _ in range(n)]
    if radians:
        angles = [numpy.rad2deg(x) for x in angles]
    return angles


def random_vectors(n=1, min_val=-100, max_val=100, vector_length=3):
    return [[(random.random() * random.choice(range(abs(max_val - min_val)))) + min_val
             for _ in range(vector_length)] for _ in range(n)]


def random_integer_vectors(n=1, min_val=-100, max_val=100, vector_length=3):
    return [[random.choice(range(min_val, max_val)) for _ in range(vector_length)] for _ in range(n)]


def random_floats(n=1, min_val=-100, max_val=100):
    return [(random.random() * random.choice(range(abs(max_val - min_val)))) + min_val for _ in range(n)]


def random_helical_helices(n=1, min_residues=10, max_residues=20):
    # build random HelicalHelix objects
    major_radii = random_floats(n=n, min_val=1, max_val=100)
    major_pitches = random_floats(n=n, min_val=100, max_val=1000)
    aas = [random.choice(range(min_residues, max_residues)) for _ in range(n)]
    phi_c_alphas = random_angles(n=n, min_val=-179, max_val=180)
    orientations = [random.choice([-1, 1]) for _ in range(n)]
    # minor_repeat can't be set to zero - must be set to None.
    minor_repeats = random_floats(n=n, min_val=0, max_val=100)
    zero_indices = [i for i, x in enumerate(minor_repeats) if x == 0.0]
    for i in zero_indices:
        minor_repeats[i] = None
    minor_helix_types = [random.choice(['alpha', 'pi', 'PPII', 'collagen']) for _ in range(n)]
    major_handedness = [random.choice(['l', 'r']) for _ in range(n)]
    hhs = [isambard.ampal.specifications.polymer_specs.helix.HelicalHelix(aa=aas[i], major_pitch=major_pitches[i],
                        major_radius=major_radii[i], major_handedness=major_handedness[i],
                        minor_helix_type=minor_helix_types[i], orientation=orientations[i],
                        phi_c_alpha=phi_c_alphas[i], minor_repeat=minor_repeats[i])
           for i in range(n)]
    return hhs

__author__ = 'Jack W. Heal'


