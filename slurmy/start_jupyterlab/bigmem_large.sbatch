#!/bin/bash
#SBATCH --time=4:0:0
#SBATCH --partition=bigmem2
#SBATCH --job-name=jlab
#SBATCH --output=nb-log-%J.out
#SBATCH --error=nb-log-%J.err
#SBATCH --mem=250G

ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)
login=midway2.rcc.uchicago.edu

# These will be written to `nb-log-%J.out`
echo -e "\n\n   Copy & paste the following to your local computer's terminal to tunnel
to the remote server:\n"
echo        "   \`\`\`bash
"
echo        "   ssh -N -L $ipnport:$ipnip:$ipnport $USER@$login -v
"
echo        "   \`\`\`
"
echo -e "\n\n   Then use a web browser on your local computer to open the following
address:\n "
echo        "   \`\`\`google-chrome
"
echo        "   localhost:$ipnport
"
echo -e     "   \`\`\`\n\n"
echo        "   You may be asked to provide a token (a string of 48 letters + digits).
Please find it in the \`nb-log-*.err\` file under your current directory.
"
sleep 1
# Start jupyter
cd /project2/gilad/
module load python cuda
jupyter-lab --no-browser --port=$ipnport --ip=$ipnip

# Then follow instructions above.
