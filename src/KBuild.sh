TARGET=$1
echo "Target chosen: " $TARGET

#Pay attention floating point array sizes
for i in 2K 4K 8K 16K 32K 64K 128K 256K 512K 5M 15M 32M
    do
      arr=$(echo $i | tr -dc '0-9.0-9')
      xplier=$(echo $i | tr -dc 'a-zA-Z')

      if [[ "${xplier}" == "K" ]]; then
          arrSz=$(echo "scale=0; ($arr * 1024)/1" | bc)
      else
          arrSz=$(echo "scale=0; ($arr * 1024 * 1024)/1" | bc)
      fi
      echo "arr=$arr, arrSz=$arrSz"

      if [ "$TARGET" = "x86_64" ] || [ "$TARGET" == "all" ]; then
          #Build x86_64 bin
          #IvyBridge

          make -f Makefile clean
          make VERBOSE=1 -f Makefile ARCH=x86_64 ArraySize=${arrSz}
          mv QStreamv1.elf bin/QStreamv1.${i}x8B.x86_64.gcc.6.3.0.linux.elf
          mv QStreamv2.elf bin/QStreamv2.${i}x8B.x86_64.gcc.6.3.0.linux.elf
          echo "Completed building x86_64 elf"

          #IvyBridge
          #make -f Makefile clean
          #make -f Makefile ARCH=x86_64 CPU=ivybridge
          #mv QStreamv1.elf bin/QStreamv1.x86_64.ivybridge.gcc.6.3.0.linux.elf
          #mv QStreamv2.elf bin/QStreamv2.x86_64.ivybridge.gcc.6.3.0.linux.elf
          #echo "Completed building x86_64 ivybridge elf"

          #Haswell
          #make -f Makefile clean
          #make -f Makefile ARCH=x86_64 CPU=haswell
          #mv QStreamv1.elf bin/QStreamv1.x86_64.haswell.gcc.6.3.0.linux.elf
          #mv QStreamv2.elf bin/QStreamv2.x86_64.haswell.gcc.6.3.0.linux.elf
          #echo "Completed building x86_64 haswell elf"

          ##Broadwell
          #make -f Makefile clean
          #make -f Makefile ARCH=x86_64 CPU=broadwell
          #mv QStreamv1.elf bin/QStreamv1.x86_64.broadwell.gcc.6.3.0.linux.elf
          #mv QStreamv2.elf bin/QStreamv2.x86_64.broadwell.gcc.6.3.0.linux.elf
          #echo "Completed building x86_64 broadwell elf"
      fi

      if [ "$TARGET" = "aarch64" ] || [ "$TARGET" == "all" ]; then
          #Build aarch64 bin
          make -f Makefile clean
          make VERBOSE=1 -f Makefile ARCH=aarch64 COMPILER=llvm ArraySize=${arrSz}
          mv QStreamv1.elf bin/QStreamv1.${i}x8B.aarch64.gcc.llvm6.0.0.linux.elf
          mv QStreamv2.elf bin/QStreamv2.${i}x8B.aarch64.gcc.llvm6.0.0.linux.elf
          echo "Completed LLVM build for" ${i}" array elements"
          echo "-----------------------------------------------"

          make -f Makefile clean
          make VERBOSE=1 -f Makefile ARCH=aarch64 COMPILER=linaro ArraySize=${arrSz}
          mv QStreamv1.elf bin/QStreamv1.${i}x8B.aarch64.gcc.lin7.1.1.linux.elf
          mv QStreamv2.elf bin/QStreamv2.${i}x8B.aarch64.gcc.lin7.1.1.linux.elf
          echo "Completed Linaro build for" ${i}" array elements"
          echo "-----------------------------------------------"

          echo "Completed building aarch64 elf"
      fi
done

