BRAF -> MAP2K1 in BioPax (NCI-PID)
==================================

* From: MEK1-2()
* Catalyzed by: BRAF(O-phospho-L-threonine@373,residue modification, active)
* To: MEK1-2(residue modification, active)
* References: http://identifiers.org/pubmed/7935374

::

    BRAF(p373) is/can be active
    BRAF(p373, active) activates Mek1 AND Mek2.

----

* From: MEK1()
* Catalyzed by: BRAF(O-phospho-L-serine@578,residue modification, active) | KSR(residue modification, active)
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References: http://identifiers.org/pubmed/19541618

.. Comment?

* From: MEK1()
* Catalyzed by: IQGAP1() | BRAF(O-phospho-L-serine@578,residue modification, active)
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References:

  * http://identifiers.org/pubmed/17563371
  * http://identifiers.org/pubmed/14970219
  * http://identifiers.org/pubmed/18567582
  * http://identifiers.org/pubmed/7731720
  * http://identifiers.org/pubmed/16135787

::

    BRAF(p578) is/can be active.
    MEK1(p222, p218) is/canbe active.
    BRAF(p578) catalyzes the reaction MEK1(???) -> MEK1(p218, p222, active)

----

* From: MEK1()
* Catalyzed by: RAF1()/BRAF()
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References: http://identifiers.org/pubmed/16508002

::

    RAF1() and BRAF() can form a complex.
    MEK1(p222, p218) is/canbe active.
    The complex RAF1()/BRAF() can catalyze the reaction
        MEK1(???) -> MEK1(p218, p222, active)

----

* From: MEK1-2()
* Catalyzed by: RAF1-BRAF(residue modification, active)
* To: MEK1-2-active(residue modification, active)
* References: http://identifiers.org/pubmed/15342917

::

    A degenerate statement? Presumably means:
    "RAF1 or BRAF catalyze activation of MEK1 or MEK2.

* From: MEK1-2()
* Catalyzed by: BRAF(residue modification, active)
* To: MEK1-2-active(residue modification, active)
* References:

  * http://identifiers.org/pubmed/17563371
  * http://identifiers.org/pubmed/7731720
  * http://identifiers.org/pubmed/16135787
  * http://identifiers.org/pubmed/7706312

::

    Included within the above.

----

* From: MEK1()
* Catalyzed by: BRAF(residue modification, active)
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References: http://identifiers.org/pubmed/10976102

::

    MEK1(p218, p222) is/canbe active.
    BRAF(active) catalyzes phosphorylation of MEK1(???) to MEK1(p218, p222, active)

----

* From: MEK1()
* Catalyzed by: BRAF(O-phospho-L-threonine@373,residue modification, active)
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References: http://identifiers.org/pubmed/7935374

.. Comment

* From: MEK1()
* Catalyzed by: BRAF(O-phospho-L-threonine@373,residue modification, active)
* To: MEK1(O-phospho-L-serine@222,O-phospho-L-serine@218,residue modification, active)
* References: http://identifiers.org/pubmed/9560161

::

    MEK1(p218, p222) is/canbe active.
    BRAF(p373) is/canbe active.
    BRAF(p373, active) catalyzes phosphorylation of MEK1(???)
        to MEK1(p218, p222, active)

----

* From: KSR(O-phospho-L-serine@297)/MEK1-2(residue modification, inactive)/14-3-3 family()
* Catalyzed by: BRAF(O-phospho-L-serine@446,O-phospho-L-threonine@598,O-phospho-L-serine@601)/14-3-3 family()
* To: KSR(O-phospho-L-serine@297)/MEK1-2-active(residue modification, active)/14-3-3 family()
* References:

  * http://identifiers.org/pubmed/7731720
  * http://identifiers.org/pubmed/7706312
  * http://identifiers.org/pubmed/19541618
  * http://identifiers.org/pubmed/7760835

::

    KSR(p297), MEK1 or MEK2(inactive), and 14-3-3 family proteins form a complex.
    BRAF(p446, p598, p601) forms a complex with 14-3-3 family proteins.
    BRAF(...)/14-3-3 catalyzes the activation of MEK1 or MEK2 in this complex.

----

* From: KSR(O-phospho-L-serine@297)/MEK1-2(residue modification, inactive)/14-3-3 family()
* Catalyzed by: 14-3-3()/BRAF(O-phospho-L-serine@446,O-phospho-L-threonine@599,O-phospho-L-serine@602)/RAF1(O-phospho-L-serine@338,O-phospho-L-serine@621,O-phospho-L-serine@494,O-phospho-L-threonine@491)
* To: KSR(O-phospho-L-serine@297)/MEK1-2-active(residue modification, active)/14-3-3 family()
* References:

  * http://identifiers.org/pubmed/19933846
  * http://identifiers.org/pubmed/19727074
  * http://identifiers.org/pubmed/7706312
  * http://identifiers.org/pubmed/19541618
  * http://identifiers.org/pubmed/12932319

::

    KSR(p297), MEK1 or MEK2(inactive), and 14-3-3 family proteins form a complex.
    BRAF(p446, p599, p602) forms a complex with RAF1(p338, p621, p494, p491)
        and 14-3-3 family proteins.
    BRAF(...)/RAF1(...)/14-3-3 catalyzes the activation of MEK1 or MEK2 in
        this complex.

Collected Factoids
------------------

::

    BRAF(p373) is/can be active
    BRAF(p373) is/canbe active.
    BRAF(p578) is/can be active.

    RAF1 or BRAF catalyze activation of MEK1 or MEK2.

    BRAF(p373, active) activates Mek1 AND Mek2.
    BRAF(p373, active) catalyzes phosphorylation of MEK1(???)
        to MEK1(p218, p222, active)
    BRAF(p578) catalyzes the reaction MEK1(???) -> MEK1(p218, p222, active)
    BRAF(active) catalyzes the reaction MEK1(???) to MEK1(p218, p222, active)
    The complex RAF1()/BRAF() can catalyze the reaction
        MEK1(???) -> MEK1(p218, p222, active)

    MEK1(p222, p218) is/canbe active.
    MEK1(p222, p218) is/canbe active.
    MEK1(p218, p222) is/canbe active.
    MEK1(p218, p222) is/canbe active.

    RAF1() and BRAF() can form a complex.
    BRAF(p446, p598, p601) forms a complex with 14-3-3 family proteins.
    KSR(p297), MEK1 or MEK2(inactive), and 14-3-3 family proteins form a complex.
    KSR(p297), MEK1 or MEK2(inactive), and 14-3-3 family proteins form a complex.

    BRAF(...)/14-3-3 catalyzes the activation of MEK1 or MEK2 in this complex.
    BRAF(p446, p599, p602) forms a complex with RAF1(p338, p621, p494, p491)
        and 14-3-3 family proteins.
    BRAF(...)/RAF1(...)/14-3-3 catalyzes the activation of MEK1 or MEK2 in
        this complex.

