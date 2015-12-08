.. _ras_minimal:

Ras (minimal)
=============

Model description
-----------------
This passage was processed and assembled automatically by INDRA:
    
    SOS1, bound to GRB2 binds RAS. RAS, bound to SOS1, binds GTP.

INDRA-assembled model component
-------------------------------

::

    def ras_minimal():
        Monomer('GTP', ['RAS'])
        Monomer(u'GRB2', [u'SOS1'])
        Monomer('RAS', [u'SOS1', 'GTP'])
        Monomer(u'SOS1', [u'GRB2', 'RAS'])

        Parameter(u'kf_sr_bind_1', 1e-06)
        Parameter(u'kr_sr_bind_1', 1e-06)
        Parameter('kf_rg_bind_1', 1e-06)
        Parameter('kr_rg_bind_1', 1e-06)
        Parameter('GTP_0', 100.0)
        Parameter(u'GRB2_0', 100.0)
        Parameter('RAS_0', 100.0)
        Parameter(u'SOS1_0', 100.0)

        Rule(u'SOS1_GRB2_RAS_bind', SOS1(GRB2=2, RAS=None) % GRB2(SOS1=2) +
            RAS(SOS1=None) <> SOS1(GRB2=2, RAS=1) % GRB2(SOS1=2) % RAS(SOS1=1),
            kf_sr_bind_1, kr_sr_bind_1)
        Rule(u'RAS_SOS1_GTP_bind', RAS(SOS1=2, GTP=None) % SOS1(RAS=2) +
            GTP(RAS=None) <> RAS(SOS1=2, GTP=1) % SOS1(RAS=2) % GTP(RAS=1),
            kf_rg_bind_1, kr_rg_bind_1)


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
