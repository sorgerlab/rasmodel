Ras regulation by RasGAPs
=========================

GTP hydrolysis by Ras is slow but is accelerated by RasGAPs.

    [9247124]_: The GTP-binding proteins return to the inactive state by virtue
    of the GTPase reaction, which is usually very slow but can be accelerated
    by the action of GAPs, in the case of the Ras/Ras-GAPs and Ran/Ran-GAP
    interactions by several orders of magnitude [1569940]_ [8262937]_
    [7548002]_.

There are several GAPs for the Ras subfamily.

    [9247124]_: Five mammalian GAPs for Ras have been described. The first,
    p120GAP, is the prototype of this class of proteins and was the first one
    to be isolated 20, 21 and 22. Apart from being a regulator of Ras, its
    N-terminal domain contains a number of signalling modules such as SH2, SH3,
    PH, Calb domains and is believed to be a signal transduction molecule that
    may act independently of Ras 23 and 24.

RASA1
-----

The domain information for RASA1 was taken from `Uniprot (ID: P20936)
<http://www.uniprot.org/uniprot/P20936>`_.

::

    #  p120GAP = RASA1.
    # CaLB (calcium lipid binding domain) is also known as a C2 domain.
    Monomer('RASA1', ['SH2_1', 'SH3', 'SH2_2', 'PH', 'C2', 'rasgap'])

NF1
---

    [9247124]_: The second Ras-specific GAP is neurofibromin (NF1), which is
    the product of the neurofibromatosis gene [25] and has also been shown to
    stimulate the GTPase of Ras 26, 27 and 28.  This gene has been found to be
    frequently mutated in patients with the disease neurofibromatosis type I
    29, 30 and 31 but also, albeit less frequently, in solid tumors [32].

`Uniprot (ID: P21359) <http://www.uniprot.org/uniprot/P21359>`_.

::

    Monomer('NF1', ['rasgap', 'CRALTRIO'])

RASA2
-----

    [9247124]_: GAP1m, a mammalian homologue of the Drosophila GAP1 gene [33],
    has been described, and a close homologue GAPIII [34], both of which
    contain, in addition to the GRD (GAP-related domain), C2 domains and a PH
    domain.

`Uniprot (ID: Q15283) <http://www.uniprot.org/uniprot/Q15283>`_.

::

    # GAP1m = RASA2
    Monomer('RASA2', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])

RASA3
-----

NOTE: As it turned out GAPIII and GAP1IP4BP are the same protein.

    [9247124]_: Recently an inositol-4-phosphate (IP4) binding protein GAP1IP4BP
    has been purified, cloned, and found to contain a Ras-GAP catalytic domain.
    In contrast to the other GAP mentioned, which are specific for Ras,
    GAP1IP4BP stimulates the GTPase of both Ras and Rap [35].

`Uniprot (ID: Q14644) <http://www.uniprot.org/uniprot/Q14644>`_.

::

    # GAPIII = GAP1IP4BP = RASA3
    Monomer('RASA3', ['C2_1', 'C2_2', 'rasgap', 'PH', 'BTK'])


