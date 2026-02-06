<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using GJP Stabilisation and dealiasing</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_GJP.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.000621718</value>
            <value variable="v" tolerance="1e-6">0.000140076</value>
            <value variable="w" tolerance="1e-6">0.000127</value>
            <value variable="p" tolerance="1e-6">0.000419658</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00202151</value>
            <value variable="v" tolerance="1e-6">0.000605264</value>
            <value variable="w" tolerance="1e-6">0.000446312</value>
            <value variable="p" tolerance="1e-6">0.00223119</value>
        </metric>
    </metrics>
</test>
