.. _raf_minimal:

Raf (minimal)
=============

Model description
-----------------
This passage was processed and assembled automatically by INDRA:

    RAF binds to the RAS-GTP complex. The RAF-RAS complex binds another RAF-RAS complex.

INDRA-assembled model component
-------------------------------

::

    def raf_minimal():
        Monomer('GTP', [u'HRAS'])
        Monomer(u'HRAS', ['GTP', u'RAF1'])
        Monomer(u'RAF1', [u'HRAS', u'RAF1'])

        Parameter(u'kf_rh_bind_1', 1e-06)
        Parameter(u'kr_rh_bind_1', 1e-06)
        Parameter(u'kf_rr_bind_1', 1e-06)
        Parameter(u'kr_rr_bind_1', 1e-06)
        Parameter('GTP_0', 100.0)
        Parameter(u'HRAS_0', 100.0)
        Parameter(u'RAF1_0', 100.0)


        Rule(u'RAF1_HRAS_GTP_bind', RAF1(HRAS=None) + HRAS(GTP=2, RAF1=None) %\
            GTP(HRAS=2) <> RAF1(HRAS=1) % HRAS(GTP=2, RAF1=1) % GTP(HRAS=2),
            kf_rh_bind_1, kr_rh_bind_1)
        Rule(u'RAF1_HRAS_RAF1_HRAS_bind', RAF1(HRAS=2, RAF1=None) %\
            HRAS(RAF1=2) + RAF1(HRAS=3, RAF1=None) % HRAS(RAF1=3) <>
            RAF1(HRAS=2, RAF1=1) % HRAS(RAF1=2) % RAF1(HRAS=3, RAF1=1) %\
            HRAS(RAF1=3), kf_rr_bind_1, kr_rr_bind_1)

.. raw:: html

    <script>
        window.setTimeout(function() {
        $('div.highlight-python pre > span.c:last-child').each(
            function () {
                if ($(this).text() == '#') {
                    $(this.nextSibling).detach();
                    $(this).detach();
                }
            }
        );
        }, 1000);
    </script>
