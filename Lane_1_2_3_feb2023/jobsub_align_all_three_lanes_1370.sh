#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=168:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_align_all_three_lanes_1370_2
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jjobsub_align_all_three_lanes_1370_2_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=48

files="R01-CR-CE-01
R01-CR-CE-02
R01-CR-CE-03
R01-CR-CE-04
R01-CR-CE-05
R01-CR-CE-06
R01-CR-CE-07
R01-CR-CE-08
R01-CR-CE-09
R01-CR-CE-10
R01-CR-CE-11
R01-CR-CE-12
R01-CR-CE-13
R01-CR-CE-14
R01-CR-CE-15
R01-CR-CE-16
R01-CR-CE-17
R01-CR-CE-18
R01-CR-CE-19
R01-CR-CE-20
R01-CR-CE-21
R01-CR-CE-22
R01-CR-CE-23
R01-CR-CE-48
R01-CR-CE-25
R01-CR-CE-26
R01-CR-CE-27
R01-CR-CE-28
R01-CR-CE-29
R01-CR-CE-30
R01-CR-RA-01
R01-CR-RA-02
R01-CR-RA-03
R01-CR-RA-04
R01-CR-RA-05
R01-CR-RA-06
R01-CR-RA-07
R01-CR-RA-08
R01-CR-RA-09
R01-CR-RA-10
R01-CR-RA-11
R01-CR-RA-12
R01-CR-RA-13
R01-CR-RA-14
R01-CR-RA-15
R01-CR-RA-16
R01-CR-RA-17
R01-CR-RA-18
R01-CR-RA-19
R01-CR-RA-20
R01-CR-RA-21
R01-CR-RA-22
R01-CR-RA-23
R01-CR-RA-48
R01-CR-RA-25
R01-CR-RA-26
R01-CR-RA-27
R01-CR-RA-28
R01-CR-RA-29
R01-CR-RA-30
R01-NI-CX-01
R01-NI-CX-02
R01-NI-CX-03
R01-NI-CX-04
R01-NI-CX-05
R01-NI-CX-06
R01-NI-CX-07
R01-NI-CX-08
R01-NI-CX-09
R01-NI-CX-10
R01-NI-CX-11
R01-NI-CX-12
R01-NI-CX-13
R01-NI-CX-14
R01-NI-CX-15
R01-NI-CX-16
R01-NI-CX-17
R01-NI-CX-18
R01-NI-CX-19
R01-NI-CX-20
R01-NI-CX-21
R01-NI-CX-22
R01-NI-CX-23
R01-NI-CX-48
R01-NI-CX-25
R01-NI-CX-26
R01-NI-CX-27
R01-NI-CX-28
R01-NI-CX-29
R01-NI-CX-30
R01-NI-OC-01
R01-NI-OC-02
R01-NI-OC-03
R01-NI-OC-04
R01-NI-OC-05
R01-NI-OC-06
R01-NI-OC-07
R01-NI-OC-08
R01-NI-OC-09
R01-NI-OC-10
R01-NI-OC-11
R01-NI-OC-12
R01-NI-OC-13
R01-NI-OC-14
R01-NI-OC-15
R01-NI-OC-16
R01-NI-OC-17
R01-NI-OC-18
R01-NI-OC-19
R01-NI-OC-20
R01-NI-OC-21
R01-NI-OC-22
R01-NI-OC-23
R01-NI-OC-48
R01-NI-OC-25
R01-NI-OC-26
R01-NI-OC-27
R01-NI-OC-28
R01-NI-OC-29
R01-NI-OC-30
R01-SI-CA-01
R01-SI-CA-02
R01-SI-CA-03
R01-SI-CA-04
R01-SI-CA-05
R01-SI-CA-06
R01-SI-CA-07
R01-SI-CA-08
R01-SI-CA-09
R01-SI-CA-10
R01-SI-CA-11
R01-SI-CA-12
R01-SI-CA-13
R01-SI-CA-14
R01-SI-CA-15
R01-SI-CA-16
R01-SI-CA-17
R01-SI-CA-18
R01-SI-CA-19
R01-SI-CA-20
R01-SI-CA-21
R01-SI-CA-22
R01-SI-CA-23
R01-SI-CA-48
R01-SI-CA-25
R01-SI-CA-26
R01-SI-CA-27
R01-SI-CA-28
R01-SI-CA-29
R01-SI-CA-30
R01-SI-DB-01
R01-SI-DB-02
R01-SI-DB-03
R01-SI-DB-04
R01-SI-DB-05
R01-SI-DB-06
R01-SI-DB-07
R01-SI-DB-08
R01-SI-DB-09
R01-SI-DB-10
R01-SI-DB-11
R01-SI-DB-12
R01-SI-DB-13
R01-SI-DB-14
R01-SI-DB-15
R01-SI-DB-16
R01-SI-DB-17
R01-SI-DB-18
R01-SI-DB-19
R01-SI-DB-20
R01-SI-DB-21
R01-SI-DB-22
R01-SI-DB-23
R01-SI-DB-48
R01-SI-DB-25
R01-SI-DB-26
R01-SI-DB-27
R01-SI-DB-28
R01-SI-DB-29
R01-SI-DB-30
R01-SI-FR-01
R01-SI-FR-02
R01-SI-FR-03
R01-SI-FR-04
R01-SI-FR-05
R01-SI-FR-06
R01-SI-FR-07
R01-SI-FR-08
R01-SI-FR-09
R01-SI-FR-10
R01-SI-FR-11
R01-SI-FR-12
R01-SI-FR-13
R01-SI-FR-14
R01-SI-FR-15
R01-SI-FR-16
R01-SI-FR-17
R01-SI-FR-18
R01-SI-FR-19
R01-SI-FR-20
R01-SI-FR-21
R01-SI-FR-22
R01-SI-FR-23
R01-SI-FR-48
R01-SI-FR-25
R01-SI-FR-26
R01-SI-FR-27
R01-SI-FR-28
R01-SI-FR-29
R01-SI-FR-30
R01-SI-GL-01
R01-SI-GL-02
R01-SI-GL-03
R01-SI-GL-04
R01-SI-GL-05
R01-SI-GL-06
R01-SI-GL-07
R01-SI-GL-08
R01-SI-GL-09
R01-SI-GL-10
R01-SI-GL-11
R01-SI-GL-12
R01-SI-GL-13
R01-SI-GL-14
R01-SI-GL-15
R01-SI-GL-16
R01-SI-GL-17
R01-SI-GL-18
R01-SI-GL-19
R01-SI-GL-20
R01-SI-GL-21
R01-SI-GL-22
R01-SI-GL-23
R01-SI-GL-48
R01-SI-GL-25
R01-SI-GL-26
R01-SI-GL-27
R01-SI-GL-28
R01-SI-GL-29
R01-SI-GL-30
R01-SI-HR-01
R01-SI-HR-02
R01-SI-HR-03
R01-SI-HR-04
R01-SI-HR-05
R01-SI-HR-06
R01-SI-HR-07
R01-SI-HR-08
R01-SI-HR-09
R01-SI-HR-10
R01-SI-HR-11
R01-SI-HR-12
R01-SI-HR-13
R01-SI-HR-14
R01-SI-HR-15
R01-SI-HR-16
R01-SI-HR-17
R01-SI-HR-18
R01-SI-HR-19
R01-SI-HR-20
R01-SI-HR-21
R01-SI-HR-22
R01-SI-HR-23
R01-SI-HR-48
R01-SI-HR-25
R01-SI-HR-26
R01-SI-HR-27
R01-SI-HR-28
R01-SI-HR-29
R01-SI-HR-30
R01-SI-KE-01
R01-SI-KE-02
R01-SI-KE-03
R01-SI-KE-04
R01-SI-KE-05
R01-SI-KE-06
R01-SI-KE-07
R01-SI-KE-08
R01-SI-KE-09
R01-SI-KE-10
R01-SI-KE-11
R01-SI-KE-12
R01-SI-KE-13
R01-SI-KE-14
R01-SI-KE-15
R01-SI-KE-16
R01-SI-KE-17
R01-SI-KE-18
R01-SI-KE-19
R01-SI-KE-20
R01-SI-KE-21
R01-SI-KE-22
R01-SI-KE-23
R01-SI-KE-48
R01-SI-KE-25
R01-SI-KE-26
R01-SI-KE-27
R01-SI-KE-28
R01-SI-KE-29
R01-SI-KE-30
R01-SI-LC-01
R01-SI-LC-02
R01-SI-LC-03
R01-SI-LC-04
R01-SI-LC-05
R01-SI-LC-06
R01-SI-LC-07
R01-SI-LC-08
R01-SI-LC-09
R01-SI-LC-10
R01-SI-LC-11
R01-SI-LC-12
R01-SI-LC-13
R01-SI-LC-14
R01-SI-LC-15
R01-SI-LC-16
R01-SI-LC-17
R01-SI-LC-18
R01-SI-LC-19
R01-SI-LC-20
R01-SI-LC-21
R01-SI-LC-22
R01-SI-LC-23
R01-SI-LC-48
R01-SI-LC-25
R01-SI-LC-26
R01-SI-LC-27
R01-SI-LC-28
R01-SI-LC-29
R01-SI-LC-30
R01-SI-LL-01
R01-SI-LL-02
R01-SI-LL-03
R01-SI-LL-04
R01-SI-LL-05
R01-SI-LL-06
R01-SI-LL-07
R01-SI-LL-08
R01-SI-LL-09
R01-SI-LL-10
R01-SI-LL-11
R01-SI-LL-12
R01-SI-LL-13
R01-SI-LL-14
R01-SI-LL-15
R01-SI-LL-16
R01-SI-LL-17
R01-SI-LL-18
R01-SI-LL-19
R01-SI-LL-20
R01-SI-LL-21
R01-SI-LL-22
R01-SI-LL-23
R01-SI-LL-48
R01-SI-LL-25
R01-SI-LL-26
R01-SI-LL-27
R01-SI-LL-28
R01-SI-LL-29
R01-SI-LL-30
R01-SI-MO-01
R01-SI-MO-02
R01-SI-MO-03
R01-SI-MO-04
R01-SI-MO-05
R01-SI-MO-06
R01-SI-MO-07
R01-SI-MO-08
R01-SI-MO-09
R01-SI-MO-10
R01-SI-MO-11
R01-SI-MO-12
R01-SI-MO-13
R01-SI-MO-14
R01-SI-MO-15
R01-SI-MO-16
R01-SI-MO-17
R01-SI-MO-18
R01-SI-MO-19
R01-SI-MO-20
R01-SI-MO-21
R01-SI-MO-22
R01-SI-MO-23
R01-SI-MO-48
R01-SI-MO-25
R01-SI-MO-26
R01-SI-MO-27
R01-SI-MO-28
R01-SI-MO-29
R01-SI-MO-30
R01-SI-SJ-01
R01-SI-SJ-02
R01-SI-SJ-03
R01-SI-SJ-04
R01-SI-SJ-05
R01-SI-SJ-06
R01-SI-SJ-07
R01-SI-SJ-08
R01-SI-SJ-09
R01-SI-SJ-10
R01-SI-SJ-11
R01-SI-SJ-12
R01-SI-SJ-13
R01-SI-SJ-14
R01-SI-SJ-15
R01-SI-SJ-16
R01-SI-SJ-17
R01-SI-SJ-18
R01-SI-SJ-19
R01-SI-SJ-20
R01-SI-SJ-21
R01-SI-SJ-22
R01-SI-SJ-23
R01-SI-SJ-48
R01-SI-SJ-25
R01-SI-SJ-26
R01-SI-SJ-27
R01-SI-SJ-28
R01-SI-SJ-29
R01-SI-SJ-30
R02-CW-AG-01
R02-CW-AG-02
R02-CW-AG-03
R02-CW-AG-04
R02-CW-AG-05
R02-CW-AG-06
R02-CW-AG-07
R02-CW-AG-08
R02-CW-AG-09
R02-CW-AG-10
R02-CW-AG-11
R02-CW-AG-12
R02-CW-AG-13
R02-CW-AG-14
R02-CW-AG-15
R02-CW-AG-16
R02-CW-AG-17
R02-CW-AG-18
R02-CW-AG-19
R02-CW-AG-20
R02-CW-AG-21
R02-CW-AG-22
R02-CW-AG-23
R02-CW-AG-48
R02-CW-AG-25
R02-CW-AG-26
R02-CW-AG-27
R02-CW-AG-28
R02-CW-AG-29
R02-CW-AG-30
R02-CW-AL-01
R02-CW-AL-02
R02-CW-AL-03
R02-CW-AL-04
R02-CW-AL-05
R02-CW-AL-06
R02-CW-AL-07
R02-CW-AL-08
R02-CW-AL-09
R02-CW-AL-10
R02-CW-AL-11
R02-CW-AL-12
R02-CW-AL-13
R02-CW-AL-14
R02-CW-AL-15
R02-CW-AL-16
R02-CW-AL-17
R02-CW-AL-18
R02-CW-AL-19
R02-CW-AL-20
R02-CW-AL-21
R02-CW-AL-22
R02-CW-AL-23
R02-CW-AL-48
R02-CW-AL-25
R02-CW-AL-26
R02-CW-AL-27
R02-CW-AL-28
R02-CW-AL-29
R02-CW-BE-01
R02-CW-BE-02
R02-CW-BE-03
R02-CW-BE-04
R02-CW-BE-05
R02-CW-BE-06
R02-CW-BE-07
R02-CW-BE-08
R02-CW-BE-09
R02-CW-BE-10
R02-CW-BE-11
R02-CW-BE-12
R02-CW-BE-13
R02-CW-BE-14
R02-CW-BE-15
R02-CW-BE-16
R02-CW-BE-17
R02-CW-BE-18
R02-CW-BE-19
R02-CW-BE-20
R02-CW-BE-21
R02-CW-BE-22
R02-CW-BE-23
R02-CW-BE-48
R02-CW-BE-25
R02-CW-BE-26
R02-CW-BE-27
R02-CW-BE-28
R02-CW-BE-29
R02-CW-BE-30
R02-CW-DE-01
R02-CW-DE-02
R02-CW-DE-03
R02-CW-DE-04
R02-CW-DE-05
R02-CW-DE-06
R02-CW-DE-07
R02-CW-DE-08
R02-CW-DE-09
R02-CW-DE-10
R02-CW-DE-11
R02-CW-DE-12
R02-CW-DE-13
R02-CW-DE-14
R02-CW-DE-15
R02-CW-DE-16
R02-CW-DE-17
R02-CW-DE-18
R02-CW-DE-19
R02-CW-DE-20
R02-CW-DE-21
R02-CW-DE-22
R02-CW-DE-23
R02-CW-DE-48
R02-CW-DE-25
R02-CW-DE-26
R02-CW-DE-27
R02-CW-DE-28
R02-CW-DE-29
R02-CW-DE-30
R02-CW-JO-01
R02-CW-JO-02
R02-CW-JO-03
R02-CW-JO-04
R02-CW-JO-05
R02-CW-JO-06
R02-CW-JO-07
R02-CW-JO-08
R02-CW-JO-09
R02-CW-JO-10
R02-CW-JO-11
R02-CW-JO-12
R02-CW-JO-13
R02-CW-JO-14
R02-CW-JO-15
R02-CW-JO-16
R02-CW-JO-17
R02-CW-JO-18
R02-CW-JO-19
R02-CW-JO-20
R02-CW-JO-21
R02-CW-JO-22
R02-CW-JO-23
R02-CW-JO-48
R02-CW-JO-25
R02-CW-JO-26
R02-CW-JO-27
R02-CW-JO-28
R02-CW-JO-29
R02-CW-JO-30
R02-CW-KA-01
R02-CW-KA-02
R02-CW-KA-03
R02-CW-KA-04
R02-CW-KA-05
R02-CW-KA-06
R02-CW-KA-07
R02-CW-KA-08
R02-CW-KA-09
R02-CW-KA-10
R02-CW-KA-11
R02-CW-KA-12
R02-CW-KA-13
R02-CW-KA-14
R02-CW-KA-15
R02-CW-KA-16
R02-CW-KA-17
R02-CW-KA-18
R02-CW-KA-19
R02-CW-KA-20
R02-CW-KA-21
R02-CW-KA-22
R02-CW-KA-23
R02-CW-KA-48
R02-CW-KA-25
R02-CW-KA-26
R02-CW-KA-27
R02-CW-KA-28
R02-CW-KA-29
R02-CW-KA-30
R02-CW-LA-01
R02-CW-LA-02
R02-CW-LA-03
R02-CW-LA-04
R02-CW-LA-05
R02-CW-LA-06
R02-CW-LA-07
R02-CW-LA-08
R02-CW-LA-09
R02-CW-LA-10
R02-CW-LA-11
R02-CW-LA-12
R02-CW-LA-13
R02-CW-LA-14
R02-CW-LA-15
R02-CW-LA-16
R02-CW-LA-17
R02-CW-LA-18
R02-CW-LA-19
R02-CW-LA-20
R02-CW-LA-21
R02-CW-LA-22
R02-CW-LA-23
R02-CW-LA-48
R02-CW-LA-25
R02-CW-LA-26
R02-CW-LA-27
R02-CW-LA-28
R02-CW-LA-29
R02-CW-LA-30
R02-CW-MI-01
R02-CW-MI-02
R02-CW-MI-03
R02-CW-MI-04
R02-CW-MI-05
R02-CW-MI-06
R02-CW-MI-07
R02-CW-MI-08
R02-CW-MI-09
R02-CW-MI-10
R02-CW-MI-11
R02-CW-MI-12
R02-CW-MI-13
R02-CW-MI-14
R02-CW-MI-15
R02-CW-MI-16
R02-CW-MI-17
R02-CW-MI-18
R02-CW-MI-19
R02-CW-MI-20
R02-CW-MI-21
R02-CW-MI-22
R02-CW-MI-23
R02-CW-MI-48
R02-CW-MI-25
R02-CW-MI-26
R02-CW-MI-27
R02-CW-MI-28
R02-CW-MI-29
R02-CW-MI-30
R02-CW-RY-01
R02-CW-RY-02
R02-CW-RY-03
R02-CW-RY-04
R02-CW-RY-05
R02-CW-RY-06
R02-CW-RY-07
R02-CW-RY-08
R02-CW-RY-09
R02-CW-RY-10
R02-CW-RY-11
R02-CW-RY-12
R02-CW-RY-13
R02-CW-RY-14
R02-CW-RY-15
R02-CW-RY-16
R02-CW-RY-17
R02-CW-RY-18
R02-CW-RY-19
R02-CW-RY-20
R02-CW-RY-21
R02-CW-RY-22
R02-CW-RY-23
R02-CW-RY-48
R02-CW-RY-25
R02-CW-RY-26
R02-CW-RY-27
R02-CW-RY-28
R02-CW-RY-29
R02-CW-SI-01
R02-CW-SI-02
R02-CW-SI-03
R02-CW-SI-04
R02-CW-SI-05
R02-CW-SI-06
R02-CW-SI-07
R02-CW-SI-08
R02-CW-SI-09
R02-CW-SI-10
R02-CW-SI-11
R02-CW-SI-12
R02-CW-SI-13
R02-CW-SI-14
R02-CW-SI-15
R02-CW-SI-16
R02-CW-SI-17
R02-CW-SI-18
R02-CW-SI-19
R02-CW-SI-20
R02-CW-SI-21
R02-CW-SI-22
R02-CW-SI-23
R02-CW-SI-48
R02-CW-SI-25
R02-CW-SI-26
R02-CW-SI-27
R02-CW-SI-28
R02-CW-SI-29
R02-CW-SI-30
R02-SC-HA-01
R02-SC-HA-02
R02-SC-HA-03
R02-SC-HA-04
R02-SC-HA-05
R02-SC-HA-06
R02-SC-HA-07
R02-SC-HA-08
R02-SC-HA-09
R02-SC-HA-10
R02-SC-HA-11
R02-SC-HA-12
R02-SC-HA-13
R02-SC-HA-14
R02-SC-HA-15
R02-SC-HA-16
R02-SC-HA-17
R02-SC-HA-18
R02-SC-HA-19
R02-SC-HA-20
R02-SC-HA-21
R02-SC-HA-22
R02-SC-HA-23
R02-SC-HA-48
R02-SC-HA-25
R02-SC-HA-26
R02-SC-HA-27
R02-SC-HA-28
R02-SC-HA-29
R02-SC-HA-30
R02-SC-IN-01
R02-SC-IN-02
R02-SC-IN-03
R02-SC-IN-04
R02-SC-IN-05
R02-SC-IN-06
R02-SC-IN-07
R02-SC-IN-08
R02-SC-IN-09
R02-SC-IN-10
R02-SC-IN-11
R02-SC-IN-12
R02-SC-IN-13
R02-SC-IN-14
R02-SC-IN-15
R02-SC-IN-16
R02-SC-IN-17
R02-SC-IN-18
R02-SC-IN-19
R02-SC-IN-20
R02-SC-IN-21
R02-SC-IN-22
R02-SC-IN-23
R02-SC-IN-48
R02-SC-IN-25
R02-SC-IN-26
R02-SC-IN-27
R02-SC-IN-28
R02-SC-IN-29
R02-SC-IN-30
R02-SS-CR-01
R02-SS-CR-02
R02-SS-CR-03
R02-SS-CR-04
R02-SS-CR-05
R02-SS-CR-06
R02-SS-CR-07
R02-SS-CR-08
R02-SS-CR-09
R02-SS-CR-10
R02-SS-CR-11
R02-SS-CR-12
R02-SS-CR-13
R02-SS-CR-14
R02-SS-CR-15
R02-SS-CR-16
R02-SS-CR-17
R02-SS-CR-18
R02-SS-CR-19
R02-SS-CR-20
R02-SS-CR-21
R02-SS-CR-22
R02-SS-CR-23
R02-SS-CR-48
R02-SS-CR-25
R02-SS-CR-26
R02-SS-CR-27
R02-SS-CR-28
R02-SS-CR-29
R02-SS-CR-30
R02-SS-FA-01
R02-SS-FA-02
R02-SS-FA-03
R02-SS-FA-04
R02-SS-FA-05
R02-SS-FA-06
R02-SS-FA-07
R02-SS-FA-08
R02-SS-FA-09
R02-SS-FA-10
R02-SS-FA-11
R02-SS-FA-12
R02-SS-FA-13
R02-SS-FA-14
R02-SS-FA-15
R02-SS-FA-16
R02-SS-FA-17
R02-SS-FA-18
R02-SS-FA-19
R02-SS-FA-20
R02-SS-FA-21
R02-SS-FA-22
R02-SS-FA-23
R02-SS-FA-48
R02-SS-FA-25
R02-SS-FA-26
R02-SS-FA-27
R02-SS-FA-28
R02-SS-FA-29
R02-SS-FA-30
R02-SS-KH-01
R02-SS-KH-02
R02-SS-KH-03
R02-SS-KH-04
R02-SS-KH-05
R02-SS-KH-06
R02-SS-KH-07
R02-SS-KH-08
R02-SS-KH-09
R02-SS-LC-01
R02-SS-LC-02
R02-SS-LC-03
R02-SS-LC-04
R02-SS-LC-05
R02-SS-LC-06
R02-SS-LC-07
R02-SS-LC-08
R02-SS-LC-09
R02-SS-LC-10
R02-SS-LC-11
R02-SS-LC-12
R02-SS-LC-13
R02-SS-LC-14
R02-SS-LC-15
R02-SS-LC-16
R02-SS-LC-17
R02-SS-LC-18
R02-SS-LC-19
R02-SS-LC-20
R02-SS-LC-21
R02-SS-LC-22
R02-SS-LC-23
R02-SS-LC-48
R02-SS-LC-25
R02-SS-LC-26
R02-SS-LC-27
R02-SS-LC-28
R02-SS-LC-29
R02-SS-LC-30
R02-SS-LT-01
R02-SS-LT-02
R02-SS-LT-03
R02-SS-LT-04
R02-SS-LT-05
R02-SS-LT-06
R02-SS-LT-07
R02-SS-LT-08
R02-SS-LT-09
R02-SS-LT-10
R02-SS-LT-11
R02-SS-LT-12
R02-SS-LT-13
R02-SS-LT-14
R02-SS-LT-15
R02-SS-LT-16
R02-SS-LT-17
R02-SS-LT-18
R02-SS-LT-19
R02-SS-LT-20
R02-SS-LT-21
R02-SS-LT-22
R02-SS-LT-23
R02-SS-LT-48
R02-SS-LT-25
R02-SS-LT-26
R02-SS-LT-27
R02-SS-LT-28
R02-SS-LT-29
R02-SS-LT-30
R06-GH-DL-01
R06-GH-DL-02
R06-GH-DL-03
R06-GH-DL-04
R06-GH-DL-05
R06-GH-DL-06
R06-GH-DL-07
R06-GH-DL-08
R06-GH-DL-09
R06-GH-DL-10
R06-GH-DL-11
R06-GH-DL-12
R06-GH-DL-13
R06-GH-DL-14
R06-GH-DL-15
R06-GH-DL-16
R06-GH-DL-17
R06-GH-DL-18
R06-GH-DL-19
R06-GH-DL-20
R06-GH-DL-21
R06-GH-DL-22
R06-GH-DL-23
R06-GH-DL-48
R06-GH-DL-25
R06-GH-DL-26
R06-GH-DL-27
R06-GH-DL-28
R06-GH-DL-29
R06-GH-DL-30
R06-GH-GW-01
R06-GH-GW-02
R06-GH-GW-03
R06-GH-GW-04
R06-GH-GW-05
R06-GH-GW-06
R06-GH-GW-07
R06-GH-GW-08
R06-GH-GW-09
R06-GH-GW-10
R06-GH-GW-11
R06-GH-GW-12
R06-GH-GW-13
R06-GH-GW-14
R06-GH-GW-15
R06-GH-GW-16
R06-GH-GW-17
R06-GH-GW-18
R06-GH-GW-19
R06-GH-GW-20
R06-GH-GW-21
R06-GH-GW-22
R06-GH-GW-23
R06-GH-GW-48
R06-GH-GW-25
R06-GH-GW-26
R06-GH-GW-27
R06-GH-GW-28
R06-GH-GW-29
R06-GH-GW-30
R06-GH-LL-01
R06-GH-LL-02
R06-GH-LL-03
R06-GH-LL-04
R06-GH-LL-05
R06-GH-LL-06
R06-GH-LL-07
R06-GH-LL-08
R06-GH-LL-09
R06-GH-LL-10
R06-GH-LL-11
R06-GH-LL-12
R06-GH-LL-13
R06-GH-LL-14
R06-GH-LL-15
R06-GH-LL-16
R06-GH-LL-17
R06-GH-LL-18
R06-GH-LL-19
R06-GH-LL-20
R06-GH-LL-21
R06-GH-LL-22
R06-GH-LL-23
R06-GH-LL-48
R06-GH-LL-25
R06-GH-LL-26
R06-GH-LL-27
R06-GH-LL-28
R06-GH-LL-29
R06-GH-LL-30
R06-GH-PQ-01
R06-GH-PQ-02
R06-GH-PQ-03
R06-GH-PQ-04
R06-GH-PQ-05
R06-GH-PQ-06
R06-GH-PQ-07
R06-GH-PQ-08
R06-GH-PQ-09
R06-GH-PQ-10
R06-GH-PQ-11
R06-GH-PQ-12
R06-GH-PQ-13
R06-GH-PQ-14
R06-GH-PQ-15
R06-GH-PQ-16
R06-GH-PQ-17
R06-GH-PQ-18
R06-GH-PQ-19
R06-GH-PQ-20
R06-GH-PQ-21
R06-GH-PQ-22
R06-GH-PQ-23
R06-GH-PQ-48
R06-GH-PQ-25
R06-GH-PT-01
R06-GH-PT-02
R06-GH-PT-03
R06-GH-PT-04
R06-GH-PT-05
R06-GH-PT-06
R06-GH-PT-07
R06-GH-PT-08
R06-GH-PT-09
R06-GH-PT-10
R06-GH-PT-11
R06-GH-PT-12
R06-GH-PT-13
R06-GH-PT-14
R06-GH-PT-15
R06-GH-PT-16
R06-GH-PT-17
R06-GH-PT-18
R06-GH-PT-19
R06-GH-PT-20
R06-GH-PT-21
R06-GH-PT-22
R06-GH-PT-23
R06-GH-PT-48
R06-GH-PT-25
R06-GH-PT-26
R06-GH-PT-27
R06-GH-PT-28
R06-GH-PT-29
R06-GH-PT-30
R06-GI-CK-01
R06-GI-CK-02
R06-GI-CK-03
R06-GI-CK-04
R06-GI-CK-05
R06-GI-CK-06
R06-GI-CK-07
R06-GI-CK-08
R06-GI-CK-09
R06-GI-CK-10
R06-GI-CK-11
R06-GI-CK-12
R06-GI-CK-13
R06-GI-CK-14
R06-GI-CK-15
R06-GI-CK-16
R06-GI-CK-17
R06-GI-CK-18
R06-GI-CK-19
R06-GI-CK-20
R06-GI-CK-21
R06-GI-CK-22
R06-GI-CK-23
R06-GI-CK-48
R06-GI-CK-25
R06-GI-CK-26
R06-GI-CK-27
R06-GI-CK-28
R06-GI-CK-29
R06-GI-CK-30
R06-GI-CN-01
R06-GI-CN-02
R06-GI-CN-03
R06-GI-CN-04
R06-GI-CN-05
R06-GI-CN-06
R06-GI-CN-07
R06-GI-CN-08
R06-GI-CN-09
R06-GI-CN-10
R06-GI-CN-11
R06-GI-CN-12
R06-GI-CN-13
R06-GI-CN-14
R06-GI-CN-15
R06-GI-CN-16
R06-GI-CN-17
R06-GI-CN-18
R06-GI-CN-19
R06-GI-CN-20
R06-GI-CN-21
R06-GI-CN-22
R06-GI-CN-23
R06-GI-CN-48
R06-GI-CN-25
R06-GI-CN-26
R06-GI-CN-27
R06-GI-CN-28
R06-GI-CN-29
R06-GI-CN-30
R06-GI-EV-01
R06-GI-EV-02
R06-GI-EV-03
R06-GI-EV-04
R06-GI-EV-05
R06-GI-EV-06
R06-GI-EV-07
R06-GI-EV-08
R06-GI-EV-09
R06-GI-EV-10
R06-GI-EV-11
R06-GI-EV-12
R06-GI-EV-13
R06-GI-EV-14
R06-GI-EV-15
R06-GI-EV-16
R06-GI-EV-17
R06-GI-EV-18
R06-GI-EV-19
R06-GI-EV-20
R06-GI-EV-21
R06-GI-EV-22
R06-GI-EV-23
R06-GI-EV-48
R06-GI-EV-25
R06-GI-EV-26
R06-GI-EV-27
R06-GI-EV-28
R06-GI-EV-29
R06-GI-EV-30
R06-GI-GL-01
R06-GI-GL-02
R06-GI-GL-03
R06-GI-GL-04
R06-GI-GL-05
R06-GI-GL-06
R06-GI-GL-07
R06-GI-GL-08
R06-GI-GL-09
R06-GI-GL-10
R06-GI-GL-11
R06-GI-GL-12
R06-GI-GL-13
R06-GI-GL-14
R06-GI-GL-15
R06-GI-GL-16
R06-GI-GL-17
R06-GI-GL-18
R06-GI-GL-19
R06-GI-GL-20
R06-GI-GL-21
R06-GI-GL-22
R06-GI-GL-23
R06-GI-GL-48
R06-GI-GL-25
R06-GI-GL-26
R06-GI-GL-27
R06-GI-GL-28
R06-GI-GL-29
R06-GI-GL-30
R06-GI-LV-01
R06-GI-LV-02
R06-GI-LV-03
R06-GI-LV-04
R06-GI-LV-05
R06-GI-LV-06
R06-GI-LV-07
R06-GI-LV-08
R06-GI-LV-09
R06-GI-LV-10
R06-GI-LV-11
R06-GI-LV-12
R06-GI-LV-13
R06-GI-LV-14
R06-GI-LV-15
R06-GI-LV-16
R06-GI-LV-17
R06-GI-LV-18
R06-GI-LV-19
R06-GI-LV-20
R06-GI-LV-21
R06-GI-LV-22
R06-GI-LV-23
R06-GI-LV-48
R06-GI-LV-25
R06-GI-LV-26
R06-GI-LV-27
R06-GI-LV-28
R06-GI-LV-29
R06-GI-LV-30
R06-GI-MY-01
R06-GI-MY-02
R06-GI-MY-03
R06-GI-MY-04
R06-GI-MY-05
R06-GI-MY-06
R06-GI-MY-07
R06-GI-MY-08
R06-GI-MY-09
R06-GI-MY-10
R06-GI-MY-11
R06-GI-MY-12
R06-GI-MY-13
R06-GI-MY-14
R06-GI-MY-15
R06-GI-MY-16
R06-GI-MY-17
R06-GI-MY-18
R06-GI-MY-19
R06-GI-MY-20
R06-GI-MY-21
R06-GI-MY-22
R06-GI-MY-23
R06-GI-MY-48
R06-GI-MY-25
R06-GI-MY-26
R06-GI-MY-27
R06-GI-MY-28
R06-GI-MY-29
R06-GI-MY-30
R06-GI-RR-01
R06-GI-RR-02
R06-GI-RR-03
R06-GI-RR-04
R06-GI-RR-05
R06-GI-RR-06
R06-GI-RR-07
R06-GI-RR-08
R06-GI-RR-09
R06-GI-RR-10
R06-GI-RR-11
R06-GI-RR-12
R06-GI-RR-13
R06-GI-RR-14
R06-GI-RR-15
R06-GI-RR-16
R06-GI-RR-17
R06-GI-RR-18
R06-GI-RR-19
R06-GI-RR-20
R06-GI-RR-21
R06-GI-RR-22
R06-GI-RR-23
R06-GI-RR-48
R06-GI-RR-25
R06-GI-RR-26
R06-GI-RR-27
R06-GI-RR-28
R06-GI-RR-29
R06-GI-RR-30
R06-MI-MM-01
R06-MI-MM-02
R06-MI-MM-03
R06-MI-MM-04
R06-MI-MM-05
R06-MI-MM-06
R06-MI-MM-07
R06-MI-MM-08
R06-MI-MM-09
R06-MI-MM-10
R06-MI-MM-11
R06-MI-MM-12
R06-MI-MM-13
R06-MI-MM-14
R06-MI-MM-15
R06-MI-MM-16
R06-MI-MM-17
R06-MI-MM-18
R06-MI-MM-19
R06-MI-MM-20
R06-MI-MM-21
R06-MI-MM-22
R06-MI-MM-23
R06-MI-MM-48
R06-MI-MM-25
R06-MI-MM-26
R06-MI-MM-27
R06-MI-MM-28
R06-MI-MM-29
R06-MI-MM-30
R06-SK-FE-01
R06-SK-FE-02
R06-SK-FE-03
R06-SK-FE-04
R06-SK-FE-05
R06-SK-FE-06
R06-SK-FE-07
R06-SK-FE-08
R06-SK-FE-09
R06-SK-FE-10
R06-SK-FE-11
R06-SK-FE-12
R06-SK-FE-13
R06-SK-FE-14
R06-SK-FE-15
R06-SK-FE-16
R06-SK-FE-17
R06-SK-FE-18
R06-SK-FE-19
R06-SK-FE-20
R06-SK-FE-21
R06-SK-FE-22
R06-SK-FE-23
R06-SK-FE-48
R06-SK-FE-25
R06-SK-FE-26
R06-SK-FE-27
R06-SK-FE-28
R06-SK-FE-29
R06-SK-FE-30
R06-SK-ON-01
R06-SK-ON-02
R06-SK-ON-03
R06-SK-ON-04
R06-SK-ON-05
R06-SK-ON-06
R06-SK-ON-07
R06-SK-ON-08
R06-SK-ON-09
R06-SK-ON-10
R06-SK-ON-11
R06-SK-ON-12
R06-SK-ON-13
R06-SK-ON-14
R06-SK-ON-15
R06-SK-ON-16
R06-SK-ON-17
R06-SK-ON-18
R06-SK-ON-19
R06-SK-ON-20
R06-SK-ON-21
R06-SK-ON-22
R06-SK-ON-23
R06-SK-ON-48
R06-SK-ON-25
R06-SK-ON-26
R06-SK-ON-27
R06-SK-ON-28
R06-SK-ON-29
R06-SK-ON-30
R06-SK-ML-01
R06-SK-ML-02
R06-SK-ML-03
R06-SK-ML-04
R06-SK-ML-05
R06-SK-ML-06
R06-SK-ML-07
R06-SK-ML-08
R06-SK-ML-09
R06-SK-ML-10
R06-SK-ML-11
R06-SK-ML-12
R06-SK-ML-13
R06-SK-ML-14
R06-SK-ML-15
R06-SK-ML-16
R06-SK-ML-17
R06-SK-ML-18"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/jobsub_align_all_three_lanes_1370_2.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/all3lanessamples/${sample}.1.fq.gz $src/Demultiplexing_stacks/all3lanessamples/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_1370_2.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20_all3lanes/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




