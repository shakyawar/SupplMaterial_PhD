package com.silicolife.MyMavenProject;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections15.ListUtils;

import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import pt.uminho.ceb.biosystems.mew.core.model.components.EnvironmentalConditions;
import silicocore.api.analysis.simulation.SimulationEnv;
import silicocore.api.container.ContainerEnv;

public class KO1analysis {
	public float getKOpredEfficiency(float threshold) throws Exception {
		final String GeneExpData = "D:\\L.major_iAC560_OptFlux\\Extended_model\\Gene_knockout_exprmnt_Data_2003.csv";
		
		// getting KOrxnList1 which are deactivated from GIMME analysis
		GimmeTests object2 = new GimmeTests();
		HashMap myMap = (HashMap) object2.getConsistency(threshold);
		
		String GIMMEDeactivatedRxnList = (String) myMap.get("key2");
		GIMMEDeactivatedRxnList = GIMMEDeactivatedRxnList.substring(1, GIMMEDeactivatedRxnList.length() - 1);
		String[] stringArray = GIMMEDeactivatedRxnList.split(",");
		ArrayList<String> GIMMEKORxnlist1 = new ArrayList<String>();
		for (int i = 0; i < stringArray.length; i++) {
			String rxn = (String) stringArray[i];
			GIMMEKORxnlist1.add(rxn);
		}
		System.out.println("GIMMEKORxnlist1 === " + GIMMEKORxnlist1);
		
		
		
		// accessing KO experimental data
		String splitBy = ",";
		BufferedReader br = new BufferedReader(new FileReader(GeneExpData));
		String line = br.readLine();
		// System.out.println(line);
		int CountWrongPrediction = 0;
		int CountCorrectPrediction = 0;

		while ((line = br.readLine()) != null) {
			String[] b = line.split(splitBy);
			String modelRxnID = b[1];
			ArrayList<String> KORxnlist2 = new ArrayList<String>(Arrays.asList(modelRxnID));
			String phenotype = b[2];
		
			// CALLING MODEL and calculating biomass value
			PromModel Promastigote = new PromModel();
			ArrayList<String> AllKOrxn = (ArrayList<String>) ListUtils.union(GIMMEKORxnlist1,KORxnlist2);
			Double Biomass = Promastigote.PromModelSimulate(AllKOrxn);
			//System.out.println("Biomass: " + Biomass);
			System.out.println("Experimental data:-- " +"Reaction ID- " + KORxnlist2 + "Phenotype- " + phenotype);


			// counting right or wrong predictions based on calculated biomass
			if (phenotype.equals("L")) {
				if (Biomass == 0.0) {
					CountCorrectPrediction = CountCorrectPrediction + 1;
					System.out.println("Correct prediction1... ");
				} else {
					CountWrongPrediction = CountWrongPrediction + 1;
					System.out.println("Wrong prediction2... ");
				}
			} else {
				if (Biomass == 0.0) {
					CountWrongPrediction = CountWrongPrediction + 1;
					System.out.println("Wrong prediction3... ");
				} else {
					CountCorrectPrediction = CountCorrectPrediction + 1;
					System.out.println("Correct prediction4... ");
				}

			}

		}

		// calculating drug target prediction efficiency
		System.out.println("-------------------Final results------------------");
		System.out.println("Total number of experimental reactions for KO analysis: "+  (CountCorrectPrediction + CountWrongPrediction));
		System.out.println("Correct predtcions = " + CountCorrectPrediction + ", Wrong predtcions = " + CountWrongPrediction);
		float predEfficiency = (CountCorrectPrediction * 100 / (CountCorrectPrediction + CountWrongPrediction));
		System.out.println("Prediction Efficiency (%) = " + predEfficiency);
		System.out.println();
		return predEfficiency;
	}
}