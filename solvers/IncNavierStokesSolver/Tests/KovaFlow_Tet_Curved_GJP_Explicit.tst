<?xml version="1.0" encoding="utf-8"?>
<test>
3   <description>3D Kovasznay flow on curved Tetrahedrons with explicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Tet_Curved.xml -I GJPStabilisation=Explicit -I GJPNormalVelocity=True</parameters>
    <files>
        <file description="Session File">KovaFlow_Tet_Curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.00028271</value>
            <value variable="v" tolerance="1e-8">0.000138304</value>
            <value variable="w" tolerance="1e-10">8.70526e-05</value>
            <value variable="p" tolerance="1e-8">0.000604777</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.00062029</value>
            <value variable="v" tolerance="1e-8">0.000503675</value>
            <value variable="w" tolerance="1e-8">0.000349651</value>
            <value variable="p" tolerance="1e-8">0.00500397</value>
        </metric>
    </metrics>
</test>
