<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Compressible Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. This test case uses the explicit solver. </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ChanFlow3D_infTurbExpl.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurbExpl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-7">0.0496262</value>
            <value variable="rhou" tolerance="1e-7">0.0510315</value>
            <value variable="rhov" tolerance="1e-8">0.00044049</value>
            <value variable="rhow" tolerance="1e-9">0.000513142</value>
            <value variable="E"    tolerance="1e-5">8.88794</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-7">0.0142084</value>
            <value variable="rhou" tolerance="1e-7">0.0220761</value>
            <value variable="rhov" tolerance="1e-8">0.00487243</value>
            <value variable="rhow" tolerance="1e-8">0.00275159</value>
            <value variable="E"    tolerance="1e-5">2.55629</value>
        </metric>
    </metrics>
</test>


