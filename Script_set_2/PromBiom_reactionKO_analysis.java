package com.silicolife.MyMavenProject;

import java.util.HashMap;


import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationProperties;
import silicocore.api.analysis.simulation.SimulationEnv;
import silicocore.api.container.ContainerEnv;

public class PromBiom {
	
	private static final String biomassEqn = "R_step_BIOMASS_LM3_sugNuc_AmastUpdated";
	private static final String modelFilePath = "D:\\L.major_iAC560_OptFlux\\Extended_model\\eiAC560_withGeneIDsIntegration.sbml";
	private static final String EnvCondition = "D:\\L.major_iAC560_OptFlux\\Integration\\GND-GFAT_media\\Medium_Amastigote_GND-GFAT.csv";
	
//	public Double getBiomass() throws Exception {
		
		
		public static void main(String[] args) throws Exception{
		
		
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath,"");
		SimulationEnv simEnv =  new SimulationEnv(c);// c.createSimulationEnv();
		simEnv.loadEnvironmentalConditions(EnvCondition, ",");
		simEnv.setTarget(biomassEqn);
		simEnv.setSimulationMethod(SimulationProperties.PFBA);
		//simEnv = simEnv.simulate();
		
		// put reaction constraints
//		simEnv = simEnv.addReactionConstraint("R_ENO", 0.00, 0.00);
		// get Biomass value
		Double flux = simEnv.getFluxValue(biomassEqn);
		System.out.println("biomass = " + flux);

		simEnv.printNetConvertion();
		
		System.out.println("ENO = " + simEnv.getFluxValue("R_ENO"));
		System.out.println (simEnv.getMapResult());
	}
	}
