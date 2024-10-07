<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Compressible Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. This test cases uses the implicit solver. </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ChanFlow3D_infTurbImpl.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurbImpl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-7">0.0496262</value>
            <value variable="rhou" tolerance="1e-7">0.0510306</value>
            <value variable="rhov" tolerance="1e-8">0.000247923</value>
            <value variable="rhow" tolerance="1e-8">0.000463459</value>
            <value variable="E"    tolerance="1e-5">8.88793</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-7">0.0140256</value>
            <value variable="rhou" tolerance="1e-7">0.0188351</value>
            <value variable="rhov" tolerance="1e-8">0.00117766</value>
            <value variable="rhow" tolerance="1e-8">0.00196196</value>
            <value variable="E"    tolerance="1e-5"> 2.51336</value>
        </metric>
    </metrics>
</test>


