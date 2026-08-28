#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/edesanti/FCCk4hep/k4geo/logs/slurm-%A_%a.out
#SBATCH --mem-per-cpu=2048M
#SBATCH -J SciFik4hep
#SBATCH -t 120:00:00

# -------------------- Protection Lines --------------------------
### Do not delete this line ###
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    ( bash "$BASH_SOURCE" "$@" )
    return $?
fi

# -------------------- Argument Parsing --------------------------
PREV_JOB=""
while [[ $# -gt 0 ]]; do
    case "$1" in
	--after)
	    PREV_JOB="$2"
	    shift 2
	    ;;
	*)
	    echo "Unknown argument: $1"
	    exit 1
	    ;;
    esac
done

# -------------------- USER INPUT PARAMETERS ----------------------
PARTICLE="mu-"
ENERGY="20*GeV"
COMMENT="SciFik4hep_singlemat_scan"

K4GEO_DIR="/home/edesanti/FCCk4hep/k4geo"
COMPACT="FCCee/ALLEGRO/compact/ALLEGRO_o2_v01/ALLEGRO_o2_v01.xml"
STEERING="FCCee/ALLEGRO/compact/ALLEGRO_o2_v01/scifi_optical_steer.py"

# Scan the gun's z-injection position along the fiber
Z_START=-1202
Z_END=1198
Z_STEP=600
SEED_BASE=3429

# n_events per z position, split into sub-runs of NEVENTSPERRUN each
NEVENTS=1000
NEVENTSPERRUN=1000

# -------------------------------------------------------------------

NZ=$(( (Z_END - Z_START) / Z_STEP + 1 ))
NRUNS=$(( NEVENTS / NEVENTSPERRUN ))
NTASKS=$(( NZ * NRUNS ))

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && -z "${SBATCH_ARRAY_TASK_RANGE:-}" ]]; then
    export SBATCH_ARRAY_TASK_RANGE="0-$((NTASKS-1))%$((NZ<48?NZ:48))"
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "COMMENT: ${COMMENT}"
    echo "NZ=${NZ}  NRUNS=${NRUNS}  NTASKS=${NTASKS}"
    echo "-------------------------------------"
    echo "Submit with:"
    COMMAND="  sbatch --array=0-$((NTASKS-1))%$((NZ<48?NZ:48)) $0"
    if [[ -n "$PREV_JOB" ]]; then
	COMMAND="${COMMAND} --dependency=after:${PREV_JOB}"
    fi
    echo $COMMAND
    exit 0
fi

# -------------------- decode this array task -----------------------
task=${SLURM_ARRAY_TASK_ID}
z_idx=$(( task % NZ ))
run_idx=$(( task / NZ ))
Z=$(( Z_START + z_idx * Z_STEP ))
SEED=$(( SEED_BASE + z_idx*100000 + run_idx ))

LOGDIR="${K4GEO_DIR}/logs"
[ -d "$LOGDIR" ] || mkdir -p "$LOGDIR"

OUTDIR="${K4GEO_DIR}/output_${PARTICLE}_${COMMENT}"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/scifi_pde_z${Z}_r${run_idx}.root"

echo "Running ddsim at gun z = $Z mm, run $run_idx, seed $SEED, $NEVENTSPERRUN events..."

# -------------------- environment ------------------------------
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
export LD_LIBRARY_PATH=${K4GEO_DIR}/build/lib:$LD_LIBRARY_PATH
pwd_save=$(pwd)
cd "${K4GEO_DIR}/InstallArea"
source bin/thisk4geo.sh
cd "${K4GEO_DIR}"
k4_local_repo InstallArea/
cd "${pwd_save}"

# -------------------- run ----------------------------------------
cd "${K4GEO_DIR}" || { echo "Failed to enter ${K4GEO_DIR}"; exit 1; }

ddsim --compactFile "${COMPACT}" \
      --steeringFile "${STEERING}" \
      --enableGun --gun.particle "${PARTICLE}" --gun.energy "${ENERGY}" \
      --gun.position "881.055*mm 0*mm ${Z}*mm" --gun.direction "1 0 0" \
      --numberOfEvents "${NEVENTSPERRUN}" \
      --random.seed "${SEED}" \
      --outputFile "${OUTFILE}"
