# reference for the installtion of the gsutil
https://cloud.google.com/sdk/docs/install#linux

# intallation of the google cli to use the gsutil
curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-400.0.0-linux-x86_64.tar.gz

tar xzfv google-cloud-cli-400.0.0-linux-x86.tar.gz

# Add the gcloud CLI to your path. Run the installation script from the root of the folder you extracted to using the following command:
./google-cloud-sdk/install.sh

# to make the command available restart the terminal or source the .bashrc file
source ~/.bashrc

# location of the installation of the 
~/Documents/software/google-cloud-sdk

# perform the authentication process make sure you manage to log from the browser
gcloud init

# command to download the segmentation information
gsutil -m cp -r \
  "gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice1/Replicate1/cell_boundaries" \
  .

# sample code to download multiple files
gsutil -m cp -r \
  "gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice2/Replicate1/cell_boundaries" \
  "gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice2/Replicate1/cell_by_gene_S2R1.csv" \
  "gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice2/Replicate1/cell_metadata_S2R1.csv" \
  "gs://public-datasets-vizgen-merfish/datasets/mouse_brain_map/BrainReceptorShowcase/Slice2/Replicate1/detected_transcripts_S2R1.csv" \
  .