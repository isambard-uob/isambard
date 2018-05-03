"""Module for evaluating the contact order of polypeptides. """

import ampal
import budeff


def calculate_contact_order(polypeptide):
    r"""Calculates the contact order of a polypeptide.

    Contact order is a is a measure of the number and range of contacts
    found in a protein normalised by sequence length [1]_. For proteins
    with folding pathways that can be approximated as two state, contact
    order is linearly related to :math:`\ln{K}` [1]_ [2]_. Contact order
    is calculated using the following method:

    .. math::

        CO = \frac{1}{LN}\sum\limits_{}^{N}\Delta{}Z_{i,j}

    Where *N* is the total number of contacts, *L* is the sequence length
    and :math:`\Delta{}Z_{i,j}` is the sequence distance between
    contacting residues.

    References
    ----------
    .. [1] Plaxco KW, Simons KT, Baker D (1998) Contact order, transition
       state placement and the refolding rates of single domain proteins,
       J. Mol Biol, **277**, 985-994.
    .. [2] Fersht AR (2000) Transition-state structure as a unifying basis
       in protein-folding mechanisms: Contact order, chain topology,
       stability, and the extended nucleus mechanism, Proc Natl Acad Sci,
       **97**, 1525-1529.

    Notes
    -----

    I've used 18 A for the Ca cut off distance to be very cautious about
    throwing away interactions. The distance between the amine and Ca of
    fully extended lysine is around 6.5 A, so if two lysines were
    interacting, it'd be 2*6.5 plus 2 times van der Waals radius, so around
    17 A.
    """
    if not isinstance(polypeptide, ampal.Polypeptide):
        raise ValueError(
            'Contact order can only be calculated for a Polypeptide chain.')
    # Force field is only assigned to enable use of find_intra_ampal
    budeff.assign_force_field(polypeptide, budeff.FORCE_FIELDS['bude_2016v1'])
    contacts = {(int(a.parent.id), int(b.parent.id))
                for a, b in budeff.find_intra_ampal(polypeptide, 18.0)
                if ampal.geometry.distance(a, b) <= 4.0}
    d_z_values = map(lambda x: x[1]-x[0]-1, contacts)
    try:
        contact_order = 1/(len(contacts)*len(polypeptide))*sum(d_z_values)
    except ZeroDivisionError:
        raise ValueError(
            'Number of contacts ({}) or length of polypeptide ({}) is 0'
            ' so contact order cannot be calculated.'.format(len(contacts),
                                                             len(polypeptide)))
    return contact_order


__author__ = 'Christopher W. Wood'
