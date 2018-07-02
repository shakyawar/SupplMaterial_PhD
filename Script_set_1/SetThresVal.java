package com.silicolife.MyMavenProject;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Dictionary;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import com.kenai.jffi.Array;
import com.sun.xml.xsom.impl.scd.Iterators.Map;

import jxl.read.biff.BiffException;
import pt.uminho.ceb.biosystems.mew.core.model.components.EnvironmentalConditions;
import pt.uminho.ceb.biosystems.mew.core.model.components.ReactionConstraint;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.gpr.ISteadyStateGeneReactionModel;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.GeneChangesList;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.GeneticConditions;
import silicocore.api.analysis.simulation.SimulationEnv;
import silicocore.api.container.ContainerEnv;
import silicocore.api.container.DataBaseContainerEnv;

import java.awt.List;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;

public class SetThresVal {

	public static void main(String[] args) throws Exception {
		Scanner thres = new Scanner(System.in);
		// System.out.println("Enter gene expression threshold value: ");
		// float user_input_threshold = thres.nextFloat();
		Set myset = new HashSet();
	
		
		for (float user_input_threshold = (float) 0; user_input_threshold < 0.5;) {
			String s = "";

			// create object of ReadGeneExpData class
			// ReadGeneExpData obj = new ReadGeneExpData();

			// Read gene expression data
			// ----------------------------------------------------------------------------------------------------------
			String splitBy = ",";
			BufferedReader br = new BufferedReader(
					new FileReader("D:\\L.major_iAC560_OptFlux\\Extended_model\\Gene_exp_data_2003.csv"));
			String line = br.readLine();
			// System.out.println(line);
			int TotalGeneExpValues = 0;
			HashMap genesLessThres_Dict = new HashMap();
			while ((line = br.readLine()) != null) {
				TotalGeneExpValues = TotalGeneExpValues + 1;
				String[] b = line.split(splitBy);
				String geneID = b[0];
				float expVal = Float.parseFloat(b[1]);
				// System.out.println(b[0]+b[1]);

				if (expVal < user_input_threshold) {
					// Using Dictionary
					//System.out.println("geneID + expVal  == " + geneID + " + " + expVal);
					genesLessThres_Dict.put(geneID, expVal);
				}
			}
			ArrayList<String> KnockoutGeneList = new ArrayList<String>(genesLessThres_Dict.keySet());
			br.close();
			// -----------------------------------------------------------------------------------------------------------

			// call method of ReadGeneExpData class using object
			// HashMap GeneDict_Selected = obj.GetData(user_input_threshold);
			System.out.println("Threshold: " + user_input_threshold);
			System.out.println();
			System.out.println("Total genes with expVal less than threshold: " + genesLessThres_Dict.size());
			System.out.println();

			// ArrayList<String> KnockoutGeneList = new
			// ArrayList<String>(GeneDict_Selected.keySet());
			System.out.println("List of genes below threshold: " + KnockoutGeneList);
			// System.out.println();
			// System.out.println();

			// // Calling model simulation from class PromModel
			// PromModel Promastigote = new PromModel();
			// // passing gene list to model
			// Double Biomass = Promastigote.PromModelSimulate(KnockoutGeneList,
			// new ArrayList<String>());
			// System.out.print("Biomass_Promastigote = " + Biomass);
			// System.out.println();


			GimmeTests object2 = new GimmeTests();
			HashMap myMap = (HashMap) object2.getConsistency(user_input_threshold);
			float inConsScore = Float.parseFloat((String) myMap.get("key_inconstScore"));
			
	//------------------------------- get reactions with threshold less and flux more vice versa ------------------------------------		
			// get count of GIMME RxnLessThres_LessFlux0 reaction list
			String GIMME_RxnLessThres_LessFlux0 = (String) myMap.get("key_RxnLessThres_LessFlux0");
			GIMME_RxnLessThres_LessFlux0 = GIMME_RxnLessThres_LessFlux0.substring(1, GIMME_RxnLessThres_LessFlux0.length() - 1);
			String[] stringArray1 = GIMME_RxnLessThres_LessFlux0.split(",");
			ArrayList<String> GIMME_RxnLessThres_LessFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray1.length; i++) {
				String rxn = (String) stringArray1[i];
				GIMME_RxnLessThres_LessFlux0_list.add(rxn);
			}
			
			String GIMME_RxnLessThres_EqualFlux0 = (String) myMap.get("key_RxnLessThres_EqualFlux0");
			GIMME_RxnLessThres_EqualFlux0 = GIMME_RxnLessThres_EqualFlux0.substring(1, GIMME_RxnLessThres_EqualFlux0.length() - 1);
			String[] stringArray2 = GIMME_RxnLessThres_EqualFlux0.split(",");
			ArrayList<String> GIMME_RxnLessThres_EqualFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray2.length; i++) {
				String rxn = (String) stringArray2[i];
				GIMME_RxnLessThres_EqualFlux0_list.add(rxn);
			}
			
			// get count of GIMME RxnLessThres_MoreFlux0 reaction list
			String GIMME_RxnLessThres_MoreFlux0 = (String) myMap.get("key_RxnLessThres_MoreFlux0");
			GIMME_RxnLessThres_MoreFlux0 = GIMME_RxnLessThres_MoreFlux0.substring(1, GIMME_RxnLessThres_MoreFlux0.length() - 1);
			String[] stringArray3 = GIMME_RxnLessThres_MoreFlux0.split(",");
			ArrayList<String> GIMME_RxnLessThres_MoreFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray3.length; i++) {
				String rxn = (String) stringArray3[i];
				GIMME_RxnLessThres_MoreFlux0_list.add(rxn);
			}			

			System.out.println("GIMME_RxnLessThres_LessFlux0_list" + GIMME_RxnLessThres_LessFlux0_list);
			System.out.println("GIMME_RxnLessThres_EqualFlux0_list" + GIMME_RxnLessThres_EqualFlux0_list);
			System.out.println("GIMME_RxnLessThres_MoreFlux0_list" + GIMME_RxnLessThres_MoreFlux0_list);
			System.out.println("---------");
		

			// get count of GIMME RxnMoreThres_MoreFlux0 reaction list
			String GIMME_RxnMoreThres_MoreFlux0 = (String) myMap.get("key_RxnMoreThres_MoreFlux0");
			GIMME_RxnMoreThres_MoreFlux0 = GIMME_RxnMoreThres_MoreFlux0.substring(1, GIMME_RxnMoreThres_MoreFlux0.length() - 1);
			String[] stringArray4 = GIMME_RxnMoreThres_MoreFlux0.split(",");
			ArrayList<String> GIMME_RxnMoreThres_MoreFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray4.length; i++) {
				String rxn = (String) stringArray4[i];
				GIMME_RxnMoreThres_MoreFlux0_list.add(rxn);
			}			
			
			// get count of GIMME RxnMoreThres_MoreFlux0 reaction list
			String GIMME_RxnMoreThres_EqualFlux0 = (String) myMap.get("key_RxnMoreThres_EqualFlux0");
			GIMME_RxnMoreThres_EqualFlux0 = GIMME_RxnMoreThres_EqualFlux0.substring(1, GIMME_RxnMoreThres_EqualFlux0.length() - 1);
			String[] stringArray5 = GIMME_RxnMoreThres_EqualFlux0.split(",");
			ArrayList<String> GIMME_RxnMoreThres_EqualFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray5.length; i++) {
				String rxn = (String) stringArray5[i];
				GIMME_RxnMoreThres_EqualFlux0_list.add(rxn);
			}
						
			// get count of GIMME RxnMoreThres_LessFlux0 reaction list
			String GIMME_RxnMoreThres_LessFlux0 = (String) myMap.get("key_RxnMoreThres_LessFlux0");
			GIMME_RxnMoreThres_LessFlux0 = GIMME_RxnMoreThres_LessFlux0.substring(1, GIMME_RxnMoreThres_LessFlux0.length() - 1);
			String[] stringArray6 = GIMME_RxnMoreThres_LessFlux0.split(",");
			ArrayList<String> GIMME_RxnMoreThres_LessFlux0_list = new ArrayList<String>();
			for (int i = 0; i < stringArray6.length; i++) {
				String rxn = (String) stringArray6[i];
				GIMME_RxnMoreThres_LessFlux0_list.add(rxn);
			}			
	//------------------------------------------------------------------------------------------------------------------------------------------------------------		
			

			// Rxn KO result for given rxn-lethalility data
//			KO1analysis object = new KO1analysis();
//			float KOPredEfficiency = object.getKOpredEfficiency(user_input_threshold);

			// System.out.println("threshold = "+ user_input_threshold);
			// System.out.println("Returned pred efficiency Score ="+
			// KOPredEfficiency);
			// System.out.println("Returned inconsitency Score ="+ inConsScore);

			s = s + user_input_threshold + "," + inConsScore + "," + GIMME_RxnLessThres_LessFlux0_list.size() + "," + GIMME_RxnLessThres_EqualFlux0_list.size() + "," + GIMME_RxnLessThres_MoreFlux0_list.size() + "," + GIMME_RxnMoreThres_MoreFlux0_list.size() + "," + GIMME_RxnMoreThres_EqualFlux0_list.size() + "," + GIMME_RxnMoreThres_LessFlux0_list.size();
			myset.add(";" + s);
			

			
			// get net conversion
			//PromModel ModelSimul = new PromModel();
			//Object NetConVal = ModelSimul.GetNetConv(KnockoutGeneList);
			//System.out.print("================Net conversion ==============" + NetConVal);
			
			
			// single gene deletion analysis
			
			PromModel ModelSimul = new PromModel();
			Double WTBiomass = ModelSimul.PromModelSimulate(KnockoutGeneList);
			
			//SignleGeneDeletion Res = new SignleGeneDeletion();
			//HashMap SingleGeneDelRes = Res.getSingGeneDelRes(user_input_threshold);
			//System.out.println("SingleGeneDelRes--------" + SingleGeneDelRes);

			DoubleGeneDeletion Res = new DoubleGeneDeletion();
			HashMap getDoubGeneDelRes = Res.getDoubGeneDelRes(user_input_threshold);
			System.out.println("getDoubGeneDelRes--------" + getDoubGeneDelRes);

			user_input_threshold = user_input_threshold + 1;
		}

		System.out.print("==============================" + myset);
	}

}