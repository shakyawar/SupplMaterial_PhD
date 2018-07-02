package com.silicolife.MyMavenProject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Spliterator;
import java.util.regex.Pattern;

import org.junit.Test;

import com.amazonaws.services.cloudfront.model.Method;

import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.core.model.components.EnvironmentalConditions;
import pt.uminho.ceb.biosystems.mew.core.model.components.ReactionConstraint;
import pt.uminho.ceb.biosystems.mew.core.model.converters.ContainerConverter;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.gpr.ISteadyStateGeneReactionModel;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.FluxValueMap;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationProperties;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SteadyStateSimulationResult;
import pt.uminho.ceb.biosystems.mew.omicsintegration.data.Condition;
import pt.uminho.ceb.biosystems.mew.omicsintegration.data.GeneDataMap;
import pt.uminho.ceb.biosystems.mew.omicsintegration.data.IOmicsContainer;
import pt.uminho.ceb.biosystems.mew.omicsintegration.data.ReactionDataMap;
import pt.uminho.ceb.biosystems.mew.omicsintegration.enums.OmicsDataType;
import pt.uminho.ceb.biosystems.mew.omicsintegration.integration.Gene2GeneIntegrator;
import pt.uminho.ceb.biosystems.mew.omicsintegration.io.CSVOmicsReader;
import pt.uminho.ceb.biosystems.mew.omicsintegration.metabolictasks.TasksReader;
import pt.uminho.ceb.biosystems.mew.omicsintegration.omicssimulation.configuration.GIMMEConfiguration;
import pt.uminho.ceb.biosystems.mew.omicsintegration.omicssimulation.methods.GIMME;
import pt.uminho.ceb.biosystems.mew.omicsintegration.othersProj.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.omicsintegration.transformation.TransformDataMapGeneToReac;
import pt.uminho.ceb.biosystems.mew.omicsintegration.utils.Config;
import pt.uminho.ceb.biosystems.mew.solvers.SolverType;
import silicocore.api.container.ContainerEnv;

public class GimmeTests {

	private static final String modelFilePath = "D:\\L.major_iAC560_OptFlux\\Extended_model\\eiAC560_withGeneIDsIntegration.sbml";
	private static final String biomassReaction = "R_step_BIOMASS_LM3_sugNuc_AmastUpdated";
	private static final String geneExpressionFilePath = "D:\\L.major_iAC560_OptFlux\\Extended_model\\Gene_exp_data_2003.csv";
	private static final String EnvCondition = "D:\\L.major_iAC560_OptFlux\\Integration\\GND-GFAT_media\\Medium_Amastigote_GND-GFAT.csv";

	
	// private static final String fluxomicsFilePath =
	// "/Users/Sara/Documents/Projects/PHD_P01_Reconst_Approaches/SimulMeythods/holm_fluxes_QP.csv";

	private static final int[] conditionsGeneExpressionDataIndex = { 1 };

	// private static final int[] conditionsFluxomicsDataIndex = {};

	// threshold for gene expression
	//private static final double[] quartile1GeneExp = { 4 };

	// private static final double[] quartile3GeneExp= { 11.41439423,
	// 11.27056071, 11.2483921 };

	// private static final double[] limitGlc = { -9.2, -11.7, -15.6 };

	// private static final String uptakeReaction =
	// "R_EX_b_DASH_D_DASH_glucose_LPAREN_e_RPAREN_";
//	private static final String biomass = "R_step_BIOMASS_LM3_sugNuc_PromastUpdated";

	private static Container container;
	private static ISteadyStateGeneReactionModel model;
	private static ArrayList<ReactionDataMap> reacsScores;
	// private static ArrayList<IOmicsContainer> fluxomicsData;
	// private static ArrayList<String> reacs;
	private List<SteadyStateSimulationResult> results;
	private List<Double> precisionList;
	private Set<String> set;
	private ArrayList<String> listDeactivatedRxn;
	private float inconstScore;
	HashMap myMap = new HashMap();

	@Test
	public Object getConsistency(float threshold) throws Exception {
		// Intialize the SBML reader for the eiAC560 Leishmania major model
		System.out.println("check point 1");
		ContainerEnv c = ContainerEnv.readSBMLFile(modelFilePath,"");
		
//		JSBMLReader sbmlReader = new JSBMLReader(modelFilePath, "Leishmania major", false);
		System.out.println("check point 2");
		container = c.getContainer(); //new Container(sbmlReader);
		// Set<String> met =
		// container.identifyMetabolitesIdByPattern(Pattern.compile(".*_b"));
		// container.removeMetabolites(met);
		container.setBiomassId(biomassReaction);
		model = (ISteadyStateGeneReactionModel) ContainerConverter.convert(container);

		// Initialize the Gene Expression Data Reader
		CSVOmicsReader geneExpressionDataReader = new CSVOmicsReader(new Condition(), geneExpressionFilePath,
				OmicsDataType.GENE);
		// geneExpressionDataReader.DELIMITER_INSIDE_FIELDS = ",";

		geneExpressionDataReader.USER_DELIMITER = ",";
		geneExpressionDataReader.setIdColumnIndex(0);
		geneExpressionDataReader.setValuesColumnIndex(1);
		geneExpressionDataReader.setHasHeader(true);

//		System.out.print(geneExpressionDataReader.getHeaders());
//		System.out.print(geneExpressionDataReader.getDiscreteValues());

		//System.out.print("check point 3");
		//System.out.println();

		// read all gene expression conditions
		reacsScores = new ArrayList<ReactionDataMap>();
		for (int i = 0; i < conditionsGeneExpressionDataIndex.length; i++) {
			geneExpressionDataReader.setValuesColumnIndex(conditionsGeneExpressionDataIndex[i]);
			IOmicsContainer dataContainer = geneExpressionDataReader.load();
			Gene2GeneIntegrator integrator = new Gene2GeneIntegrator(container, Config.FIELD_ID, Config.FIELD_ID);
			GeneDataMap expEvidence = (GeneDataMap) integrator.convert(dataContainer);

			// TransformDataMapGeneToReac transformDataMap = new
			// TransformDataMapGeneToReac(container);
			// E-FULX
			Map<String, Object> properties = new HashMap<String, Object>();
			properties.put(TransformDataMapGeneToReac.VAR_CONTAINER, container);
			properties.put(TransformDataMapGeneToReac.VAR_OPERATION_AND, "MIN");
			properties.put(TransformDataMapGeneToReac.VAR_OPERATION_OR, "PLUS");
			System.out.println("check point 4");
			TransformDataMapGeneToReac transformDataMap = new TransformDataMapGeneToReac(properties);
			reacsScores.add((ReactionDataMap) transformDataMap.transform(expEvidence));

		}

//		System.out.println("Reaction scores = " + reacsScores.get(0).getMapValues().size());

		System.out.println("check point 5");
		
		// Initialize the Fluxomics Data Reader
		// CSVOmicsReader fluxesDataReader = new CSVOmicsReader(new Condition(),
		// fluxomicsFilePath, OmicsDataType.REACTION);
		// fluxesDataReader.DELIMITER_INSIDE_FIELDS = ";";
		// fluxesDataReader.USER_DELIMITER = ",";
		// fluxesDataReader.setIdColumnIndex(0);
		// fluxesDataReader.setHasHeader(true);
		// fluxomicsData = new ArrayList<IOmicsContainer>();
		// for (int i = 0; i < conditionsFluxomicsDataIndex.length; i++) {
		// fluxesDataReader.setValuesColumnIndex(conditionsFluxomicsDataIndex[i]);
		// IOmicsContainer dataContainer = fluxesDataReader.load();
		// fluxomicsData.add(dataContainer);
		// }
		// reacs = new ArrayList<String>();
		// reacs.addAll(fluxomicsData.get(0).getValues().keySet());

		// Collections.sort(reacs);

		// for (String r : reacs)
		// System.out.println(r);

		results = new ArrayList<SteadyStateSimulationResult>();
		precisionList = new ArrayList<Double>();
		System.out.println("\n --------------- GIMME--------------- ");

		GIMMEConfiguration config = new GIMMEConfiguration(container, SolverType.CPLEX3);
		TasksReader tasks = new TasksReader(container, "D:\\L.major_iAC560_OptFlux\\Extended_model\\Task.txt",
				Config.FIELD_ID, true);
		System.out.println("check point 6");
		tasks.load();
		
		config.setTasks(tasks.getTasks());
		config.setRMFPercentage(0.999);

		for (int i = 0; i < conditionsGeneExpressionDataIndex.length; i++) {
			// Set Condition

			config.setCutOff(threshold);

			// setting environmental condition
			String splitBy = ",";
			BufferedReader br = new BufferedReader(new FileReader(EnvCondition));
			String line = br.readLine();
			// System.out.println(line);
			EnvironmentalConditions environmentalConditions = new EnvironmentalConditions();
			while ((line = br.readLine()) != null) {
				String[] b = line.split(splitBy);
				String rxnID = b[0];
				float lowBound = Float.parseFloat(b[1]);
				float upBound = Float.parseFloat(b[2]);
				environmentalConditions.addReactionConstraint(rxnID, new ReactionConstraint(lowBound, upBound));
				// System.out.println(b[0]+b[1]+b[1]);
			}
			br.close();
			// created by Sara
			// EnvironmentalConditions environmentalConditions = new
			// EnvironmentalConditions();
			// environmentalConditions.addReactionConstraint(uptakeReaction, new
			// ReactionConstraint(limitGlc[i], 0));
			// environmentalConditions.addReactionConstraint("R_EX_ile_DASH_L_LPAREN_e_RPAREN_",
			// new ReactionConstraint(-34, 0));

			System.out.println("Environmental conditions: " + environmentalConditions);
			config.setEnvironmentalConditions(environmentalConditions);
			config.setReactionScores(reacsScores.get(i));

			// run gimme
			GIMME gimme = new GIMME(model, config);
			gimme.setProperty(SimulationProperties.METHOD_NAME, SimulationProperties.PFBA);
			
			
			SteadyStateSimulationResult result = gimme.simulate();
			
			System.out.println(gimme.getMethod());
			
			//System.out.println("Otpmial Solution: " + result.getOFvalue());
			results.add(result);
			System.out.println("----------------");

			// Getting flux values results & calculating inconsistency score

			set = model.getReactions().keySet();
			System.out.println("list of all reactions - " + set);
			Object[] rxnArray = set.toArray();
			float inconstScore = (float) 0.0;

			
			int count_RxnLessThres_LessFlux0 = 0;
			ArrayList<String> RxnLessThres_LessFlux0 =  new ArrayList<String>();
			
			int count_RxnLessThres_EqualFlux0 = 0;
			ArrayList<String> RxnLessThres_EqualFlux0 =  new ArrayList<String>();
			
			int count_RxnLessThres_MoreFlux0 = 0;
			ArrayList<String> RxnLessThres_MoreFlux0 =  new ArrayList<String>();
			
			
			
			int count_RxnMoreThres_MoreFlux0 = 0;
			ArrayList<String> RxnMoreThres_MoreFlux0 =  new ArrayList<String>();
			
			int count_RxnMoreThres_EqualFlux0 = 0;
			ArrayList<String> RxnMoreThres_EqualFlux0 =  new ArrayList<String>();
			
			int count_RxnMoreThres_LessFlux0 = 0;
			ArrayList<String> RxnMoreThres_LessFlux0 =  new ArrayList<String>();
			
		
			for (int j = 0; j < rxnArray.length; j++) {
				String rxn = (String) rxnArray[j];
				Double fluxVal = result.getFluxValues().getValue(rxn);
			//	 System.out.println(rxn + "mapped val = " + reacsScores.get(0).getValue(rxn));
	//			System.out.println(rxn +"==="+ fluxVal);
				// mapped flux values is reacsScores.get(0).getValue(rxn))
				Double mapRxnVal = reacsScores.get(0).getValue(rxn);

				if (mapRxnVal == null) {
					// None

				} else if (Double.isNaN(mapRxnVal)) {
					// None				
				} 
				
				//-------------- for gene expression value > threshold -------------------	
				else if (mapRxnVal > threshold) {
					
					if (fluxVal > 0) {
						count_RxnMoreThres_MoreFlux0 = count_RxnMoreThres_MoreFlux0 + 1;
						RxnMoreThres_MoreFlux0.add(rxn);
					}
					
					else if (fluxVal == 0) {
							count_RxnMoreThres_EqualFlux0 = count_RxnMoreThres_EqualFlux0 + 1;
							RxnMoreThres_EqualFlux0.add(rxn);
					}
					
					else {
						count_RxnMoreThres_LessFlux0 = count_RxnMoreThres_LessFlux0 + 1;
						RxnMoreThres_LessFlux0.add(rxn);
					}
				}
				
				
				else {
				//-------------- for gene expression value < threshold -------------------
					if (mapRxnVal < threshold) {
						if (fluxVal > 0) {
							//System.out.println(rxn + " = mappedValue = " + mapRxnVal + " = FluxValue = " + fluxVal);
							inconstScore = (float) ((float) inconstScore
									+ (fluxVal * (float) ((threshold - mapRxnVal))));
							count_RxnLessThres_MoreFlux0 = count_RxnLessThres_MoreFlux0 + 1;
							RxnLessThres_MoreFlux0.add(rxn);
						} 
						
						else if (fluxVal == 0) {
							count_RxnLessThres_EqualFlux0 = count_RxnLessThres_EqualFlux0 + 1;
							RxnLessThres_EqualFlux0.add(rxn);
						}
						
						else {
							//System.out.println(rxn + " = mappedValue = " + mapRxnVal + " = FluxValue = " + fluxVal);
							count_RxnLessThres_LessFlux0 = count_RxnLessThres_LessFlux0 + 1;
							RxnLessThres_LessFlux0.add(rxn);
						}

					} else {
						// None
					}
				}
			}

			// created by Sara for getting final flux values of all reactions
			//System.out.println("------------------------------Final flux values of all reactions---------------------------");
			for (String r : model.getReactions().keySet())
				 System.out.println(r + " = " +
				 result.getFluxValues().getValue(r));

			System.out.println("-------------Consistency analysis ---------- ");
			System.out.println("Total number of reactions (with gene expression < threshold) and non-zero flux = " + Integer(count_RxnLessThres_MoreFlux0, count_RxnLessThres_LessFlux0));
			System.out.println("These are reactions that may reactivate = " + RxnLessThres_MoreFlux0 + RxnLessThres_LessFlux0);
			System.out.println("----------------------- ");

			System.out.println("Total number of reactions (with gene expression < threshold) and zero flux = " + count_RxnLessThres_EqualFlux0);
			System.out.println("These are reactions that are deactivated = " + RxnLessThres_EqualFlux0);
			System.out.println("----------------------- ");
			
			System.out.println("Inconsistency Score = " + inconstScore);
			System.out.println("Otpmial Solution (By Sara) = " + result.getOFvalue());
			System.out.println("Biomass = " + result.getFluxValues().getValue(biomassReaction));
			//return inconstScore;
			myMap.put("key_inconstScore", Float.toString(inconstScore));
			myMap.put("key_RxnLessThres_LessFlux0", RxnLessThres_LessFlux0.toString());
			myMap.put("key_RxnLessThres_MoreFlux0", RxnLessThres_MoreFlux0.toString());
			myMap.put("key_RxnMoreThres_MoreFlux0", RxnMoreThres_MoreFlux0.toString());
			myMap.put("key_RxnMoreThres_LessFlux0", RxnMoreThres_LessFlux0.toString());
			myMap.put("key_RxnMoreThres_EqualFlux0", RxnMoreThres_EqualFlux0.toString());
			myMap.put("key_RxnLessThres_EqualFlux0", RxnLessThres_EqualFlux0.toString());
			
			
			
			
			
			
		    return myMap;

		}
		return null;
	}


private int Integer(int i, int j) {
		// TODO Auto-generated method stub
	int k = i+j;	
		return k;
	}


//	public ArrayList<String> getGIMMEReactionsCounts(float threshold) throws Exception {
//		GimmeTests Obj = new GimmeTests();
//		Object countGIMMErxns = Obj.getConsistency(threshold);
//		
//	
//		return countGIMMErxns;
//	}
	
	
	public ArrayList<String> getGIMMEDeactiavtedGeneList(float threshold) throws Exception {
		String splitBy = ",";
		BufferedReader br = new BufferedReader(
				new FileReader(geneExpressionFilePath));
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
			if (expVal < threshold) {
				// Using Dictionary
				//System.out.println("geneID + expVal  == " + geneID + " + " + expVal);
				genesLessThres_Dict.put(geneID, expVal);
			}
		}
		ArrayList<String> GIMMEKnockoutGeneList = new ArrayList<String>(genesLessThres_Dict.keySet());
		br.close();
		System.out.println("GIMME KO removing genes  == " + GIMMEKnockoutGeneList);
		return GIMMEKnockoutGeneList;

		
	}
	}
	  	
 

