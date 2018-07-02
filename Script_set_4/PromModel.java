package com.silicolife.MyMavenProject;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections15.ListUtils;

import com.google.common.base.Objects;

import analysis.utils.SimulationUtils;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.StoichiometryValueCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.gpr.ISteadyStateGeneReactionModel;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.GeneChangesList;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.GeneticConditions;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.ReactionChangesList;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationProperties;
import silicocore.api.analysis.simulation.SimulationEnv;
import silicocore.api.container.ContainerEnv;
import silicocore.api.container.DataBaseContainerEnv;

public class PromModel {
	
	private static final String biomassEqn = "R_step_BIOMASS_LM3_sugNuc_AmastUpdated";
	private static final String modelFilePath = "D:\\L.major_iAC560_OptFlux\\Extended_model\\eiAC560_withGeneIDsIntegration.sbml";
	private static final String EnvCondition = "D:\\L.major_iAC560_OptFlux\\Integration\\GND-GFAT_media\\Medium_Amastigote_GND-GFAT.csv";

	
	public Double PromModelSimulate(ArrayList<String> KORxnlist2) throws Exception {
		return PromModelSimulate(new ArrayList<String>(), KORxnlist2);
	}

	public Double PromModelSimulate(Collection<String> geneList, ArrayList<String> KORxnlist2) throws Exception {
		
		// accessing model from database
		//DataBaseContainerEnv c_prom = DataBaseContainerEnv.loadDBContainerFromServer("pvilaca","!#prcv1l4c4", "iCA506", 1, "");
		//c_prom.addPathwayFromFile("D:\\L.major_iAC560_OptFlux\\Integration\\NGlycan_LPG_GIPL_reaction_6.txt", "c", "f");
		//c_prom.addPathwayFromFile("D:\\L.major_iAC560_OptFlux\\Integration\\NGlycan_LPG_GIPL_reaction_6.txt","r","NGlycan_LPG_GIPL_reaction_6",true, true, true);
		
		//System.out.print("------------------ Database connection started -----------------------------");
		//System.out.println();
		
		
		// accessing from the saved sbml file
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath,"");
		
		// getting the reactions that can be knocked out from the model based on GPR rules and given gene list
		ArrayList<String> KORxnlist1 = new ArrayList<String> ();
		if(geneList == null){
			geneList = new ArrayList<String>();
		}
		else {
			GeneChangesList geneList1 = new GeneChangesList(geneList);
			GeneticConditions geneCond = new GeneticConditions(geneList1, (ISteadyStateGeneReactionModel) c.getModel(), false);
			KORxnlist1.addAll(geneCond.getReactionList().getReactionIds());
			//ArrayList<String> KORxnlist1 = new ArrayList<String>(geneCond.getReactionList().getReactionIds());
		}
		
				
		//System.out.print("KORxnlist1: " + KORxnlist1.getClass().getSimpleName());
		//System.out.print("KORxnlist2: " + KORxnlist2.getClass().getSimpleName());

		
		List<String> AllKOrxn = ListUtils.union(KORxnlist1,KORxnlist2);
		//System.out.print("Total Knocked out reactions  "+ AllKOrxn);
		//System.out.println();
	
		if( AllKOrxn.size() < 1 ) {
	         //System.out.print("No reactions are knocked out ");
	         //System.out.println();
	      }else {
	    	  c.removeReactions(AllKOrxn);;
	    	  //System.out.print("Reactions are knocked out ");
	    	  //System.out.println();
	      }
		//System.out.println("Knocked out reactions (mapped from gene list): " + KORxnlist1);
		//System.out.println("Knocked out reactions (mapped from reaction list): "+ KORxnlist2);
		

		// Simulating model
		//System.out.println("Model simulation started...");
		SimulationEnv simEnv =  new SimulationEnv(c);// c.createSimulationEnv();
		simEnv.loadEnvironmentalConditions(EnvCondition, ",");
		simEnv.setTarget(biomassEqn);
		simEnv.setSimulationMethod(SimulationProperties.PFBA);
		simEnv = simEnv.simulate();
		// get Biomass value
		Double flux = simEnv.getFluxValue(biomassEqn);
		//System.out.println("Flux = " + flux);	
	return 	flux;	
	}
	
	
	
	// method to get All genes in the model
	public Set<String> GetGeneList() throws Exception {
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath,"");
		return c.getContainer().getGenes().keySet();
	}
	
	
	
	// get flux values for biomass precursors
	
	public HashMap GetBiomPrecurFlux(ArrayList<String> RxnList) throws Exception {
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath, "");
		c.removeReactions(RxnList);
		SimulationEnv simEnv = new SimulationEnv(c);
		simEnv.loadEnvironmentalConditions(EnvCondition, ",");
		simEnv.setTarget(biomassEqn);
		simEnv.setSimulationMethod(SimulationProperties.PFBA);
		simEnv = simEnv.simulate();
		
		// get biomass precursor
		ReactionChangesList reactionList = new ReactionChangesList(RxnList);
		GeneticConditions geneCond = new GeneticConditions(reactionList);
		Set<StoichiometryValueCI> metab = new HashSet<StoichiometryValueCI>(c.getContainer().getReaction(biomassEqn).getReactants().values());
		Map<String, Double> BiomPrecurFlux = SimulationUtils.testMetabolitesProductionByStoichiometry(c.getContainer(), metab, simEnv.getEc(), geneCond);
		return (HashMap) BiomPrecurFlux;
	}
	
	
	
	// get reactions that can be knocked out for given genes, based on GPR rules
	public ArrayList<String> GetRxnGPR(List<String> list) throws Exception {
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath,"");
		HashMap BiomPrecurFlux = new HashMap();
		GeneChangesList geneList1 = new GeneChangesList(list);
		GeneticConditions geneCond = new GeneticConditions(geneList1, (ISteadyStateGeneReactionModel) c.getModel(),
				false);
		ArrayList<String> KORxnlist2 = new ArrayList<String>(geneCond.getReactionList().getReactionIds());
		return KORxnlist2;
	}
	
	
	
	// get reactions precursors 
	public HashMap GetRxnPrecur(ArrayList<String> AllKOrxn) throws Exception {
		HashMap<String, String> BiomPrecurFlux = new HashMap<String, String>();
		//HashMap  BiomPrecurFlux = new HashMap();
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath, "");
		c.removeReactions(AllKOrxn);
		SimulationEnv simEnv = new SimulationEnv(c);
		simEnv.loadEnvironmentalConditions(EnvCondition, ",");
		simEnv.setTarget(biomassEqn);
		simEnv.setSimulationMethod(SimulationProperties.PFBA);
		simEnv = simEnv.simulate();
		
		//RxnPrecur = simEnv.testReactionPrecursors("R_PYK");
		
		return null;	
	
	}
	
	// get Net conversion 
	public SimulationEnv GetNetConv(ArrayList<String> AllKOrxn) throws Exception {
		//HashMap  BiomPrecurFlux = new HashMap();
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath, "");
		c.removeReactions(AllKOrxn);
		SimulationEnv simEnv = new SimulationEnv(c);
		simEnv.loadEnvironmentalConditions(EnvCondition, ",");
		simEnv.setTarget(biomassEqn);
		simEnv.setSimulationMethod(SimulationProperties.PFBA);
		simEnv = simEnv.simulate();
		SimulationEnv NetConv = simEnv.printNetConvertion();
		
		//RxnPrecur = simEnv.testReactionPrecursors("R_PYK");
		
		return NetConv;	
	
	}
	
	
	
	
}