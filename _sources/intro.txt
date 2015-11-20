Ras Executable Model Overview
=============================

The Ras Executable Model (REM) is a model of the signaling pathway of the Ras
family proteins (K-RAS, N-RAS, H-RAS), including their regulators and
effectors, at a biochemical level of detail. It is a `rule-based
model <http://www.nature.com/nmeth/journal/v8/n2/full/nmeth0211-130.html>`_
written in `PySB. <http://www.pysb.org>`_

The Ras Executable Model is being developed as part of `DARPA's Big Mechanism
program <http://www.darpa.mil/Our_Work/I2O/Programs/Big_Mechanism.aspx>`_ by the
`Sorger Lab <http://sorrger.med.harvard.edu>`_ at Harvard Medical School and
collaborators within the Big Mechanism program. In the context of the
program's overall effort to automate the construction of models from
the scientific literature, the (manual) development of the Ras model has the
following motivations:

1. Reference standard for evaluating reading algorithms.
2. Resource for improving reading algorithms.
3. Systems biology resource.
4. Investigation of the manual reading and modeling process.
5. Opportunity to engage biological domain experts.

A "literate" model
------------------

To better understand the relationship between scientific texts and executable
models, in REM the documentation and evidence supporting the various
biochemical reactions is a primary, rather than a supporting feature. As a
means of more fully integrating textual evidence with executable modeling code
we have adopted a `"literate" style of programming.
<http://en.wikipedia.org/wiki/Literate_Programming>`_  As first described by
computer scientist Donald Knuth, in literate programming,

    Instead of imagining that our main task is to instruct a computer what to
    do, let us concentrate rather on explaining to human beings what we want a
    computer to do.

Our goal for REM is for a human reader with biological and/or modeling
expertise to read it, much like a scientific paper or review, and evaluate
whether the implementation of the model's various mechanistic assertions are
adequately supported by the textual evidence.

To do this, we exploit the fact that in PySB, biological models are implemented
as Python programs, which allows tools and conventions for computer programming
to be applied to the organization and documentation of models. In REM, we embed
the PySB modeling code within the documentation itself, which is written in
`reStructuredText. <http://docutils.sourceforge.net/rst.html>`_ The
documentation is processed and formatted using the Python documentation
generator, `Sphinx, <http://sphinx-doc.org>`_ while the remainder of the build
process extracts the model code from the documentation, and executes it to
generate various outputs: exports of the model in several formats (BNGL, Kappa,
SBML), visualizations, and simulation results.

Ty Thomson's `Yeast Pheromone Model <http://yeastpheromonemodel.org>`_ is an
excellent prior example of this general approach.

Installation
------------

To generate the text of the documentation, `Sphinx <http://sphinx-doc.org>`_ is
the only dependency, which can be installed via Pip with::

    pip install -U Sphinx

To extract and run the PySB model, `PySB <http://pysb.org>`_ is required, along
with its various dependencies (including Numpy, Scipy, SymPy, Matplotlib,
BioNetGen, and optionally `KaSim <http://github.com/Kappa-Dev/KaSim>`_).

Team members and contacts
-------------------------

RAS Model (REM)
~~~~~~~~~~~~~~~
* John Bachman
* Benjamin Gyori
* Jeremy Muhlich
* Kartik Subramanian

Additional INDRA developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Daniel Carlin
* William Chen

Principal investigators
~~~~~~~~~~~~~~~~~~~~~~~
* Peter Sorger
* Dexter Pratt

Natural language parsers
~~~~~~~~~~~~~~~~~~~~~~~~
* TRIPS - James Allen, IHMC / University of Rochester
* REACH - Mihai Surdeanu, University of Arizona

If you have questions or would like to contribute, please contact
`John Bachman <http://github.com/johnbachman>`_.

Acknowledgments
---------------

This work was supported by the DARPA Big Mechanism Program under Contract No.
W911NF-14-1-0397 "Programmatic modelling for reasoning across complex
mechanisms," Peter Sorger, William Chen and Dexter Pratt PIs.

