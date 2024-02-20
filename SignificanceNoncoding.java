/************************************************************           
 * Noncoding Significance Analysis of Hg19					*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2021 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This script executes all subscripts of this 	*
 * 			significance analysis, prepares input files		*
 * 			for the subscripts, downloads annotation files	*
 * 			if needed and deletes intermediate files. The	*
 * 			parameters of this script and the format of		*
 * 			the input files are described in the user 		*
 * 			manual provided with this method.				*
 * 															*   
 ************************************************************/

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class SignificanceNoncoding {
	
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};//,,"Y"
	static String[] chr2={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"};//,,"Y"
	static String separator="/";
	static String out_suffix="";
	static ArrayList<java.util.regex.Pattern> disabled_tests=new ArrayList<java.util.regex.Pattern>();
	static int do_significance = 0;   // NB. <1 ~ no, 0 ~ automatic, >1 ~ yes
	static int do_combine = 0;

	static boolean delete_intermediate=false;

	private static java.util.Random rng = new java.util.Random();

	public static long[] Random_nextLongs(java.util.Random rng, int count) {
		long[] values = new long[count];
		for (int i = 0; i < count; ++i)
			values[i] = rng.nextLong();

		return values;
	}
	
	//parse arguments of the method and decide whether annotation files need to be downloaded
	public static void main(String[] args){
		try{
		delete_intermediate=true;
		ArrayList<String> arg=new ArrayList<String>();
		for (int i=0;i<args.length;i++){
			if(!args[i].equals("")){
				if(args[i].equals("-k")){
					delete_intermediate=false;
				}
				else if(args[i].startsWith("-z")){
					String level=args[i].substring("-z".length());
					if (!level.equals("")){
						ZipFilter.gzip_level=Integer.parseInt(level);
					}
					out_suffix=".gz";
				}
				else if(args[i].equals("-disable_test")){
					disabled_tests.add(globPattern(args[++i]));
				}
				else if(args[i].equals("-no_significance")){
					do_significance=-1;
				}
				else if(args[i].equals("-always_significance")){
					do_significance=+1;
				}
				else if(args[i].equals("-no_combine")){
					do_combine=-1;
				}
				else if(args[i].equals("-always_combine")){
					do_combine=+1;
				}
				else if(args[i].startsWith("-seed=")){
					rng = new java.util.Random(Long.parseLong(args[i].substring("-seed=".length())));
				}
				else{
					arg.add(args[i]);	
				}
					
			}
		}		

		if(arg.size()==3){
			execute(arg.get(0), arg.get(1), "", arg.get(2), delete_intermediate, true);
		}
		else if(arg.size()==4){
			execute(arg.get(0), arg.get(1), arg.get(3), arg.get(2), delete_intermediate, false);
		}

		}
		catch(java.lang.Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	//this method combines the execution of all subscripts, reformats input files for the subscripts and downloads annotation data if needed 
	public static void execute(String entity, String file_list_mutation_files, String folder_annotation, String folder_auxiliary, boolean delete_intermediate, boolean download_annotations) throws java.io.IOException {
		boolean isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");
		if(isWindows){
			separator="\\";
		}
		else{
			separator="/";
		}
		if(!new File(file_list_mutation_files).exists()){
			System.out.println(file_list_mutation_files+" does not exist");
			return;
		}
		
		if(!new File(folder_auxiliary).exists()){
			new File(folder_auxiliary).mkdirs();
		}
		
		if(download_annotations){
			folder_annotation=download_unzip(folder_auxiliary);
		}
		
		
		if(!ZipOverlay.exists(new File(folder_annotation))){
			System.out.println(folder_annotation+" does not exist");
			return;
		}
		
		folder_annotation=new File(folder_annotation).getAbsolutePath()+separator;
		folder_auxiliary=new File(folder_auxiliary).getAbsolutePath()+separator;
		
		String folder_significance=folder_auxiliary+"Significance"+separator;
		String folder_counts_all=folder_auxiliary+"CountsAll"+separator;
		
		if(!new File(folder_significance).exists()){
			new File(folder_significance).mkdir();
		}
		if(!new File(folder_counts_all).exists()){
			new File(folder_counts_all).mkdir();
		}

		if(!new File(folder_auxiliary+"MutationFiles"+separator).exists()){
			new File(folder_auxiliary+"MutationFiles"+separator).mkdir();
		}
		
		{
			ArrayList<String> entities=new ArrayList<String>();
			ArrayList<ArrayList<String>> files=new ArrayList<ArrayList<String>>();
			FileInputStream in=new FileInputStream(file_list_mutation_files);
			BufferedReader input= new BufferedReader(new InputStreamReader(ZipFilter.filterInputStream(in)));
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(!contains(t[0],entities)){
					entities.add(t[0]);
					files.add(new ArrayList<String>());
				}
				files.get(index(t[0],entities)).add(transform_path(t[1],file_list_mutation_files));
			}
			input.close();
			
			
			String[][][] files_mut_snv=new String[entities.size()][chr.length][1];
			String[][][] files_mut_indel=new String[entities.size()][chr.length][1];
			String[][] files_donors=new String[entities.size()][1];

			String [] all_entities=to_array(entities);

			boolean input_done = true;
			for (int i=0;i<files.size();i++){
				for (int ii=0;ii<chr.length;ii++){
					files_mut_snv[i][ii][0]=folder_auxiliary+"MutationFiles"+separator+"Mutations_Chr"+chr[ii]+"_"+entities.get(i)+"_SNV.txt"+out_suffix;
					input_done = input_done && new java.io.File(files_mut_snv[i][ii][0]).exists();

					files_mut_indel[i][ii][0]=folder_auxiliary+"MutationFiles"+separator+"Mutations_Chr"+chr[ii]+"_"+entities.get(i)+"_Indel.txt"+out_suffix;
					input_done = input_done && new java.io.File(files_mut_indel[i][ii][0]).exists();
				}

				files_donors[i][0]=folder_auxiliary+"MutationFiles"+separator+"Donors_"+entities.get(i)+".txt"+out_suffix;
				input_done = input_done && new java.io.File(files_donors[i][0]).exists();
			}

			System.out.println("A");
			if (input_done) {
				System.out.println("[<< reusing previous MutationFiles>>]");
			}
			else {
			for (int i=0;i<files.size();i++){
				System.out.println(files.get(i));
				
				BufferedWriter[][] output= new BufferedWriter[chr.length][2];
				ArrayList<String> donors=new ArrayList<String>();
				for (int ii=0;ii<chr.length;ii++){
					for (int jj=0;jj<2;jj++){
						if(jj==0){
							String file=folder_auxiliary+"MutationFiles"+separator+"Mutations_Chr"+chr[ii]+"_"+entities.get(i)+"_SNV.txt"+out_suffix;
							files_mut_snv[i][ii][0]=file;
							java.io.FileOutputStream out=new java.io.FileOutputStream(file);
							output[ii][jj]= new BufferedWriter(new java.io.OutputStreamWriter(ZipFilter.filterOutputStream(out, file)));
							
						}
						else{
							String file=folder_auxiliary+"MutationFiles"+separator+"Mutations_Chr"+chr[ii]+"_"+entities.get(i)+"_Indel.txt"+out_suffix;
							files_mut_indel[i][ii][0]=file;
							java.io.FileOutputStream out=new java.io.FileOutputStream(file);
							output[ii][jj]= new BufferedWriter(new java.io.OutputStreamWriter(ZipFilter.filterOutputStream(out, file)));
						}
					}
				}
				
				files_donors[i][0]=folder_auxiliary+"MutationFiles"+separator+"Donors_"+entities.get(i)+".txt"+out_suffix;
				java.io.FileOutputStream out_donor=new java.io.FileOutputStream(files_donors[i][0]);
				BufferedWriter output_donor= new BufferedWriter(new java.io.OutputStreamWriter(ZipFilter.filterOutputStream(out_donor, files_donors[i][0])));
		
				for (int j=0;j<files.get(i).size();j++){
					in=new FileInputStream(files.get(i).get(j));
					input= new BufferedReader(new InputStreamReader(ZipFilter.filterInputStream(in)));
					String[] header=input.readLine().split("	");
					
					int index_chr=index("Chromosome",header);
					int index_pos=index("Position",header);
					int index_ref=index("Reference_Allele",header);
					int index_alt=index("Tumor_Seq_Allele",header);
					int index_sample=index("Tumor_Sample_Barcode",header);
					
					
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						int ii=index(t[index_chr],chr,chr2);
						int jj=-1;
						if(t[index_ref].length()==1&&t[index_alt].length()==1){
							jj=0;
						}
						else{
							jj=1;
						}
						if(ii!=-1){
							output[ii][jj].write(t[index_sample]+"	"+t[index_pos]+"	"+t[index_ref]+"	"+t[index_alt]);
							output[ii][jj].newLine();
							if(!contains(t[index_sample],donors)){
								donors.add(t[index_sample]);
							}
						}
					
					}
					input.close();
					
				}
				
				for (int j=0;j<donors.size();j++){
					output_donor.write(donors.get(j));
					output_donor.newLine();
				}
				
				for (int ii=0;ii<output.length;ii++){
					for (int jj=0;jj<output[ii].length;jj++){
						output[ii][jj].close();
					}
				}
				output_donor.close();
			}
			} // !input_done
			System.out.println("A");
			
			{
				long seeds[] = Random_nextLongs(rng, 4);
			for (int i=0;i<4;i++){
				if (do_significance == 0 ? !new File(folder_auxiliary + separator + "Significance" + separator + String.format("Significance_%s_%d.txt", entity, i*2500) + out_suffix).exists() : do_significance > 0){
					System.out.println(i);
					CombinedStatistics_10.execute(entity, i*2500, folder_auxiliary, folder_significance, folder_annotation, folder_counts_all,  all_entities, files_donors,  files_mut_snv,  files_mut_indel, new java.util.Random(seeds[i]));
				}
			}
			}
			{
				long seeds[] = Random_nextLongs(rng, 4);
			for (int i=0;i<4;i++){
				if (do_significance == 0 ? !new File(folder_auxiliary + separator + "Significance" + separator + String.format("Significance_100_%s_%d.txt", entity, i*25000) + out_suffix).exists() : do_significance > 0){
					System.out.println(i);
					CombinedStatistics_100.execute(entity, i*25000, folder_auxiliary,  folder_significance, folder_annotation, folder_counts_all,  all_entities, files_donors, files_mut_snv,  files_mut_indel, new java.util.Random(seeds[i]));
				}
			}
			}
			
			{
				long seed = rng.nextLong();
			if (do_combine == 0 ? !new File(folder_auxiliary + separator + String.format("FDR_Weighted_Combined_%s.txt", entity) + out_suffix).exists() : do_combine > 0){
				Combine_PValues_FDR.execute(entity, folder_annotation, folder_significance,  folder_auxiliary, new java.util.Random(seed));
			}
			}
		
			if(delete_intermediate){
				if(download_annotations){
					delete_annotation(folder_annotation);
				}		
				delete(folder_auxiliary+"Significance"+separator);
				delete(folder_auxiliary+"CountsAll"+separator);
				delete(folder_auxiliary+"MutationFiles"+separator);
					
				for (int i=0;i<4;i++){
					delete(folder_auxiliary+entity+"_10kb_"+i*2500+separator);
					delete(folder_auxiliary+entity+"_100kb_"+i*25000+separator);
				}
			

			}
		
		}
	}
	
	//if the path file lists only file names and not absolute paths, this method adds the full path
	public static String transform_path(String file, String folder){
		if(new File(file).exists()){
			return file;
		}
		else if (new File(new File(folder).getParent()+separator+file).exists()){
			return new File(folder).getParent()+separator+file;
		}
		else{
			System.out.println(file+" not found");
			System.exit(0);
			return "";
		}
		
	}
	
	//delete the annotation folder
	public static void delete_annotation(String folder){
		delete(folder+"Align36mer_Dichotomous"+separator);
		delete(folder+"ASAnnotationHg19"+separator);
		delete(folder+"CoverageFiles"+separator);
		delete(folder+"ExcludeRegions"+separator);
		delete(folder+"Expression"+separator);
		delete(folder+"GeneTrackIGV_New"+separator);
		delete(folder+"Hg19"+separator);
		delete(folder+"SNV_Raw_10"+separator);
		delete(folder+"SummarizedSignal"+separator);
		new File(folder).delete();
	}
	
	
	//delete intermediate files
	public static void delete(String folder){
		System.out.println("delete "+folder);
		File[] files=new File(folder).listFiles();
		for (int i=0;i<files.length;i++){
			new File(files[i].getAbsolutePath()).delete();
		}
		new File(folder).delete();
	}
	
	//download and unzip annotation folder
	public static String download_unzip(String folder_destination) throws java.io.IOException {
		folder_destination=new File(folder_destination).getAbsolutePath()+separator;
			BufferedInputStream inputStream = new BufferedInputStream(new URL("http://storage.googleapis.com/noncoding_analysishg19/AnnotationFilesComplete.zip").openStream());
			FileOutputStream fileOS = new FileOutputStream(folder_destination+"AnnotationFilesComplete.zip");
			byte data[] = new byte[1024];
			int byteContent;
			while ((byteContent = inputStream.read(data, 0, 1024)) != -1) {
				fileOS.write(data, 0, byteContent);
			}
		
		new File(folder_destination+"AnnotationFilesComplete"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"Align36mer_Dichotomous"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"ASAnnotationHg19"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"CoverageFiles"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"ExcludeRegions"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"Expression"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"GeneTrackIGV_New"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"Hg19"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"SNV_Raw_10"+separator).mkdir();
		new File(folder_destination+"AnnotationFilesComplete"+separator+"SummarizedSignal"+separator).mkdir();
		
		FileInputStream fis;
        byte[] buffer = new byte[1024];
        {
            fis = new FileInputStream(folder_destination+"AnnotationFilesComplete.zip");
            ZipInputStream zis = new ZipInputStream(fis);
            ZipEntry ze = zis.getNextEntry();
            while(ze != null){
                String fileName = ze.getName();
                File newFile = new File(new File(folder_destination).getAbsolutePath()+separator+fileName);
        
                if(newFile.isDirectory()){
                	if(!newFile.exists()){
                		newFile.mkdir();
                	}
                }
                else{
                	FileOutputStream fos = new FileOutputStream(newFile);
                    int len;
                    while ((len = zis.read(buffer)) > 0) {
                    	fos.write(buffer, 0, len);
                    }
                    fos.close();
                }
  
                zis.closeEntry();
                ze = zis.getNextEntry();
            }
            zis.closeEntry();
            zis.close();
            fis.close();
            
            new File(folder_destination+"AnnotationFilesComplete.zip").delete();
            
        }
	        
	    return folder_destination+"AnnotationFilesComplete"+separator;
		
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static String[] to_array(ArrayList<String> t){
		String[] tt=new String[t.size()];
		for (int i=0;i<t.size();i++){
			tt[i]=t.get(i);
		}
		return tt;
	}
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	
	public static int index(String s, String[] t1, String[] t2){
		for (int i=0;i<t1.length;i++){
			if(t1[i].equals(s)){
				return i;
			}
		}
		for (int i=0;i<t2.length;i++){
			if(t2[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	public static int index (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}

	private static java.util.regex.Pattern globPattern(String pattern){
		// NB. I guess not perfect but matches Pattern.quote()
		return java.util.regex.Pattern.compile("^\\Q" + pattern.replace("*", "\\E.*\\Q").replace("?", "\\E.\\Q") + "\\E$");
	}
}
