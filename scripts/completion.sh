#/usr/bin/env bash

_sbewh_completions()
{
    local additionalKeys=("-h" "-help" "-par")
    cword="${COMP_WORDS[$COMP_CWORD]}"
    helpMode=0
    if [ $COMP_CWORD -gt 0 ]; then
        prev=${COMP_WORDS[$COMP_CWORD-1]}
        if ( [ "${prev}" = "-help" ] || [ "${prev}" = "-h" ]) ; then
            helpMode=1
            if  [ "${cword:0:1}" != "-" ] ; then
                cword="-${cword}"
            fi
        fi
    fi
    if [ "${cword:0:1}" != "-" ] ; then
        if [ ${helpMode} -eq 0 ] ; then
            COMPREPLY=( $(compgen -f -- "$cword"; compgen -d -S / -- "$cword") )
        fi
        return
    fi
    group=$(echo ${cword} | awk -F '.' '{print $1}')
    if [ ${#group} -eq ${#cword} ] ; then
        COMPREPLY=( $(compgen -W "-Diagonalization -FourierTransform -General -KspaceRegions -Output -Parallelization -Propagation -Pulse -TightBinding -h -help -par" -- $cword) )
    fi
    if [ ${#COMPREPLY[@]} -eq 1 ] ; then
        if [[ ${additionalKeys[@]} =~ $cword ]] ; then
            return
        fi
        if [ ${helpMode} -eq 0 ] ; then
            cword="${COMPREPLY[0]}."
        fi
    fi
    group=$(echo ${cword} | awk -F '.' '{print $1}')
      case ${group} in
        -Diagonalization)
            COMPREPLY=( $(compgen -W "-Diagonalization.maxRelError -Diagonalization.maxSweepCount -Diagonalization.mode" -- $cword) );;
        -FourierTransform)
            COMPREPLY=( $(compgen -W "-FourierTransform.Nf -FourierTransform.fftPlanningMode -FourierTransform.maxNfPrime -FourierTransform.performDirectCalculation -FourierTransform.planningTimeLimit -FourierTransform.saveIndices -FourierTransform.whTransformTestField" -- $cword) );;
        -General)
            COMPREPLY=( $(compgen -W "-General.dimensionality -General.dryRun -General.finishMultiRuns -General.runMode -General.timeIt -General.verbose -General.veryVerbose" -- $cword) );;
        -KspaceRegions)
            COMPREPLY=( $(compgen -W "-KspaceRegions.densityOutputStride -KspaceRegions.storeDensityMatrix" -- $cword) );;
        -Output)
            COMPREPLY=( $(compgen -W "-Output.disableRunSaving -Output.outputDirectory -Output.run -Output.saveAsNpz -Output.saveInput -Output.saveWannierInput" -- $cword) );;
        -Parallelization)
            COMPREPLY=( $(compgen -W "-Parallelization.independentMpiIntegrations -Parallelization.maxThreadCount -Parallelization.rmaPollTime -Parallelization.useGPU" -- $cword) );;
        -Propagation)
            COMPREPLY=( $(compgen -W "-Propagation.Nk -Propagation.T2 -Propagation.allowNkAdaption -Propagation.epsAbs -Propagation.epsRel -Propagation.kGlobalShift -Propagation.occupationSmearingWidth -Propagation.outputCount -Propagation.soothingWidth" -- $cword) );;
        -Pulse)
            COMPREPLY=( $(compgen -W "-Pulse.Emax -Pulse.Emax0 -Pulse.Emax1 -Pulse.Emax2 -Pulse.Emax3 -Pulse.Emax4 -Pulse.Emax5 -Pulse.Emax6 -Pulse.Emax7 -Pulse.Emax8 -Pulse.Emax9 -Pulse.FWHM -Pulse.FWHM0 -Pulse.FWHM1 -Pulse.FWHM2 -Pulse.FWHM3 -Pulse.FWHM4 -Pulse.FWHM5 -Pulse.FWHM6 -Pulse.FWHM7 -Pulse.FWHM8 -Pulse.FWHM9 -Pulse.cep -Pulse.cep0 -Pulse.cep1 -Pulse.cep2 -Pulse.cep3 -Pulse.cep4 -Pulse.cep5 -Pulse.cep6 -Pulse.cep7 -Pulse.cep8 -Pulse.cep9 -Pulse.cyclesOn -Pulse.cyclesOn0 -Pulse.cyclesOn1 -Pulse.cyclesOn2 -Pulse.cyclesOn3 -Pulse.cyclesOn4 -Pulse.cyclesOn5 -Pulse.cyclesOn6 -Pulse.cyclesOn7 -Pulse.cyclesOn8 -Pulse.cyclesOn9 -Pulse.lambda -Pulse.lambda0 -Pulse.lambda1 -Pulse.lambda2 -Pulse.lambda3 -Pulse.lambda4 -Pulse.lambda5 -Pulse.lambda6 -Pulse.lambda7 -Pulse.lambda8 -Pulse.lambda9 -Pulse.multiRunFile -Pulse.pol -Pulse.pol0 -Pulse.pol1 -Pulse.pol2 -Pulse.pol3 -Pulse.pol4 -Pulse.pol5 -Pulse.pol6 -Pulse.pol7 -Pulse.pol8 -Pulse.pol9 -Pulse.risingCycles -Pulse.risingCycles0 -Pulse.risingCycles1 -Pulse.risingCycles2 -Pulse.risingCycles3 -Pulse.risingCycles4 -Pulse.risingCycles5 -Pulse.risingCycles6 -Pulse.risingCycles7 -Pulse.risingCycles8 -Pulse.risingCycles9 -Pulse.tCentral -Pulse.tCentral0 -Pulse.tCentral1 -Pulse.tCentral2 -Pulse.tCentral3 -Pulse.tCentral4 -Pulse.tCentral5 -Pulse.tCentral6 -Pulse.tCentral7 -Pulse.tCentral8 -Pulse.tCentral9 -Pulse.tEnd -Pulse.tEnd0 -Pulse.tEnd1 -Pulse.tEnd2 -Pulse.tEnd3 -Pulse.tEnd4 -Pulse.tEnd5 -Pulse.tEnd6 -Pulse.tEnd7 -Pulse.tEnd8 -Pulse.tEnd9 -Pulse.tStart -Pulse.tStart0 -Pulse.tStart1 -Pulse.tStart2 -Pulse.tStart3 -Pulse.tStart4 -Pulse.tStart5 -Pulse.tStart6 -Pulse.tStart7 -Pulse.tStart8 -Pulse.tStart9 -Pulse.type -Pulse.type0 -Pulse.type1 -Pulse.type2 -Pulse.type3 -Pulse.type4 -Pulse.type5 -Pulse.type6 -Pulse.type7 -Pulse.type8 -Pulse.type9" -- $cword) );;
        -TightBinding)
            COMPREPLY=( $(compgen -W "-TightBinding.fermiLevel -TightBinding.occupiedBands -TightBinding.occupiedBelowBandGap -TightBinding.onlyRealMatrixElements -TightBinding.symmetrizeHamiltonian -TightBinding.temperature -TightBinding.wannierSeed" -- $cword) );;
  esac

    for key in "${additionalKeys[@]}" ; do
        if [[ ${COMP_WORDS[@]} =~ "${key}" ]] ; then
            return
        fi
    done
    if [ ${#COMPREPLY[@]} -eq 1 ] && [ ${helpMode} -eq 0 ] ; then
        COMPREPLY=( "${COMPREPLY[0]}=")
    fi
}

complete -F _sbewh_completions -o nospace ./sbewh
