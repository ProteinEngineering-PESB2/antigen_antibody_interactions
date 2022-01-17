import pandas as pd
from Bio import SeqIO, SeqRecord, Seq
import subprocess, os, time
from pandas import json_normalize
from multiprocessing import Pool
import sys
def get_values_GO(filename, dictionary, go_type, go_abb):
    file_now=open(filename, "r")
    all_lines=file_now.read().split("\n")
    if len(all_lines) > 1:
        first_result= all_lines[0].split("\t")
        dictionary[go_type+"_GO"]= first_result[1]
        dictionary["Predict_value_"+go_abb]= first_result[2]
    else:
	    dictionary[go_type+"_GO"]= "Not results"
	    dictionary["Predict_value_"+go_abb]= "-"
    return dictionary

def process_go(dict_all, path_out):
    command= "metastudent -i {} -o {}seq_GO_".format(path_out, path_out.replace("fasta/", "temp/").replace("fasta",""))
    salida=subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    files=os.listdir("./temp")
    files = [value for value in files if "_GO_" in value]
    sub_dict = {}
    get_values_GO("./temp/"+files[0], sub_dict, "Biological_Process", "BPO")
    get_values_GO("./temp/"+files[1], sub_dict, "Celular_Component", "CCO")
    get_values_GO("./temp/"+files[2], sub_dict, "Molecular_Function", "MFO")
    os.system("rm ./temp/{} ./temp/{} ./temp/{}".format(files[0], files[1], files[2]))
    dict_all["Gene Ontology"] = sub_dict
    return dict_all

def process_pfam(dict_all, fasta):
    command= "curl -k -LH 'Expect:' -F seq='<{}' -F output=xml 'https://pfam.xfam.org/search/sequence'".format(fasta)					
    salida=subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    link_search=str(salida).split("result_url")
    array_pfam_responses = []
    if len(link_search) > 1:
        link_result=link_search[1].replace(">","").replace("<","")
        link_result=link_result[:-1]
        link_completo= "https://pfam.xfam.org"+link_result
        #print(link_completo)
        time.sleep(30)
        command2= "curl -k -s -LH 'Expect:' '"+link_completo+"'"
        result = os.popen(command2).read()
        #print(result)
        pfam_file_name = "{}_pfam_result.xml".format(fasta.replace("fasta/", "temp/"))
        file_pfam= open(pfam_file_name, "w")
        file_pfam.write(result)
        file_pfam.close()
        file_pfam_in= open(pfam_file_name, "r")
        all_lines=file_pfam_in.read().split("\n")
        #print(all_lines)
        if '    <matches>' in all_lines:
            for n in range(len(all_lines)):
                if "match accession" in all_lines[n]:
                    dic_response_pfma = {"Accession": "", "Id_accession":"", "Type": "", "Class":"", "Evalue":"", "Bitscore":""}
                    #row_out=[dataset_in["Id_sequence"][i]]
                    linea1= all_lines[n].replace("          ","").replace(">","").replace("<","").split(" ")
                    linea2= all_lines[n+1].replace("            ","").replace(">","").replace("<","").split(" ")
                    dic_response_pfma["Accession"]=(linea1[1].replace('"',"").split("=")[1])
                    dic_response_pfma["Id_accession"]=(linea1[2].replace('"',"").split("=")[1])
                    dic_response_pfma["Type"]=(linea1[3].replace('"',"").split("=")[1])
                    dic_response_pfma["Class"]=(linea1[4].replace('"',"").split("=")[1])
                    dic_response_pfma["Evalue"]=(linea2[7].replace('"',"").split("=")[1])
                    dic_response_pfma["Bitscore"]=(linea2[8].replace('"',"").split("=")[1])
                    array_pfam_responses.append(dic_response_pfma)
            dict_all.update({"pfam_predicts":array_pfam_responses})
        else:
            dict_all.update({"pfam_predicts":array_pfam_responses})            
    else:
        dict_all.update({"pfam_predicts":array_pfam_responses})
    try:
        os.system("rm " + pfam_file_name)
    except:
        pass
    return dict_all

def process_structural(dict_all, fasta):
    ss3_properties=["H_ss3", "E_ss3", "C_ss3"]
    ss8_properties= [ "H_ss8", "G_ss8", "I_ss8", "E_ss8", "B_ss8", "T_ss8", "S_ss8", "L_ss8"]
    acc_properties=["B_acc", "M_acc", "E_acc"]
    disso_properties= [".", "*"]
    sh_file_properties = "Predict_Property/Predict_Property.sh"
    
    #create command
    path_out = fasta.replace("fasta/", "temp/") + "_structural"
    command= "{} -i {} -o {}".format(sh_file_properties, fasta, path_out)
    os.system(command)
    #process all results
    files=os.listdir(path_out)
    first = next(filter(lambda files: ".all" in files, files), None)
    #lecture file
    file_now=open(path_out+"/"+first, "r")
    all_lines=file_now.read().split("\n")
    file_now.close()
    sub_dict = {}
    #get properties count
    sub_dict.update({"counts_ss3":get_properties(ss3_properties, all_lines[2]), "counts_ss8": get_properties(ss8_properties, all_lines[3]), "counts_acc": get_properties(acc_properties, all_lines[4]), "counts_disorder":get_properties(disso_properties, all_lines[5])})

    #get results
    sub_dict.update({"ss3_properties":all_lines[2], "diso_properties":all_lines[5], "ss8_properties":all_lines[3], "sac_properties":all_lines[4]})
    try:
        os.system("rm -r {} {}".format(path_out, fasta))
    except:
        pass
    dict_all["Structural"] = sub_dict
    return dict_all
    #export data	
    #dict_response_for_sequence.update({"structural_predict": dict_all})

def get_properties(array_prop, array_seq):
	dict_response = {value: array_seq.count(value.split("_")[0]) for value in array_prop}
	return dict_response

def process(row):
    row = row[1]
    id = row[0]
    print("procesando", row[0])
    try:
        record = SeqRecord.SeqRecord(Seq.Seq(row.sequence), row.id, "", "")
        path_out = "fasta/{}_temp.fasta".format(id)
        SeqIO.write(record, path_out, "fasta")
        dict_all={"id": row.id, "sequence": row.sequence}
        dict_all = process_go(dict_all, path_out)
        dict_all = process_pfam(dict_all, path_out)
        dict_all = process_structural(dict_all, path_out)
        f = open("Success.txt", "a")
        f.write(row[0] + "\n")
        f.close()
	
        print("exitoso", row[0])
    except:
        f = open("Errors.txt", "a")
        f.write(row[0] + "\n")
        f.close()
        dict_all = {"id": row.id, "sequence": row.sequence}
        print("fallido", row[0])
    return dict_all
if __name__ == '__main__':
    try:
        os.mkdir("temp")
    except:
        pass
    try:
        os.mkdir("fasta")
    except:
        pass
    file = sys.argv[1]
    data = pd.read_csv(file)
    characterized = []
    with Pool(4) as p:
        characterized = p.map(process, data.iterrows())
    os.system("rm -r fasta temp")
    df = json_normalize(characterized)
    df.to_csv("{}_result.csv".format(file.replace(".csv", "")), index=False)
