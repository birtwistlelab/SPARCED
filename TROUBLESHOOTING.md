## Troubleshooting

### "Repository corrupted" or "Newer revision available"
 `kube-runner/kube-login.sh <pvc-name>` - the workflow has been updated since your last run, so manually delete the workflow contents stored in the PVC (which Nextflow stores in the `projects` folder)

### Docker `var` directory issues
A known problem with Docker on Mac is running Docker containers that make any use of the `/var/` directory (there's a symlink involved, look it up at your own peril). To fix this and run the pipeline locally on Mac, you have to change your default Docker settings. Go to Docker->Preferences->File Sharing and remove `/private` from your shared directories. Then add `/private/var/folders` and `/var/folders`. This should fix any related issues.

### SPARCED-nf gives unexpected output
Verify you have the proper input data -> `kube-runner/kube-login.sh <pvc-name>`

 - Using `kube-login.sh <pvc-name` to obtain a shell into the PVC is also helpful for checking workflow output in any of the `.command` files nextflow generates in its `workDir` during runtime
