import ROOT
import os

# List of root files to be merged
root_files = [
     "v22r.root",
     "v2all.root",
    "v22b.root",
    "v2n.root",
    "v24r.root",
    "v24b.root",
    "v24p.root",
    "exprougev2.root",
    "massjpsi.root",
    "masspsi2s.root",
    "massr.root",
    "massb.root",
    "massp.root",

]

# Output file name
output_file = "merged_outputexp.root"

# Open the output file
output = ROOT.TFile(output_file, "RECREATE")

# Loop through all the root files and merge them
for file_name in root_files:
    input_file = ROOT.TFile.Open(file_name)
    input_file.cd()
    base_name = os.path.splitext(os.path.basename(file_name))[0]

    for key in input_file.GetListOfKeys():
        obj = key.ReadObj()
        new_name = f"{key.GetName()}_{base_name}"

        if obj.IsA().InheritsFrom("TTree"):
            if not output.GetListOfKeys().Contains(new_name):
                output.cd()
                obj_new = obj.CloneTree(-1)
                obj_new.SetName(new_name)
                obj_new.Write()
            else:
                output.cd()
                obj_new = output.Get(new_name)
                obj_new.Merge(obj)
                obj_new.Write()
        else:
            if not output.GetListOfKeys().Contains(new_name):
                output.cd()
                obj_new = obj.Clone()
                obj_new.SetName(new_name)
                obj_new.Write()
            else:
                output.cd()
                obj_new = output.Get(new_name)
                if obj.IsA().InheritsFrom("TH1"):
                    obj_new.Add(obj)
                obj_new.Write()

    input_file.Close()

# Close the output file
output.Close()

print(f"Merged {len(root_files)} files into {output_file}")
