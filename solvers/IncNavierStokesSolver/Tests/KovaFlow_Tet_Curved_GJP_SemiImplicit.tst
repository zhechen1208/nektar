<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Kovasznay flow on curved Tetrahedrons with Semi-Implicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Tet_Curved.xml -I GJPStabilisation=SemiImplicit</parameters>
    <files>
        <file description="Session File">KovaFlow_Tet_Curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-14">0.000644841</value>
            <value variable="v" tolerance="1e-14">0.000234722</value>
            <value variable="w" tolerance="1e-14">0.000131382</value>
            <value variable="p" tolerance="1e-14">0.00203655</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-14">0.00505111</value>
            <value variable="v" tolerance="1e-14">0.0027092</value>
            <value variable="w" tolerance="1e-14">0.000683692</value>
            <value variable="p" tolerance="1e-14">0.0102901</value>
        </metric>
    </metrics>
</test>
