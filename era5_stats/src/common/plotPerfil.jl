#=
Project BRAMS-FURNAS
Author: Rodrigues, L.F. [LFR]
e-mail: luizfrodrigues@protonmail.com

This program read a CSV output of BRAMS' model, and extract the wind,
and than plot a vertical data in animated form to use in Furnas-BRAMs page

You will need the Julia Language installed, get in https://julialang.org/downloads/
and You must install the packages. To install the packages
use the command using Pkg and Pkg.add("Package Name") inside Julia.

Packages: CSV, DataFrames, Plots, ArgParse, ProgressBars

License: CC-GPL 3.0

How to use:

julia plotPerfil.jl <name_of_csv>

ex:

julia plotPerfil.jl Pico_Pirineus_v_20230605-TI240h-IT0000-DT30.0s-FR30.0s.csv

=#
using CSV
using DataFrames
using Plots
using ArgParse
using ProgressBars

@time begin
   file = ARGS[1]
   println("File input: ",file)

   totlines = 0

   println("Adjusting header of CSV...")
   lines = readlines(file) #Le o arquivo CSV
   f = open(file*".work","w")
   iline = 0
   #Percorre as linhas do arquivo para ajustar o cabeçalho
   for line in lines
      global iline = iline+1
      #Descarta as primeiras duas linhas
      if iline<3
         continue
      end
      #Na terceira linha escreve um cabeçalho e 
      #uma linha com as alturas
      if iline == 3
         write(f,line*"\n")
         #Substitui o m das alturas por vazios (vira número)
         write(f,replace(line, "m" => "")*"\n")
         continue
      end
      global totlines = totlines+1
      # Escreve no arquivo de trabalho
      write(f,line*"\n")
   end
   close(f)
   println("Number of times: ",totlines)

   # Abre o arquivo de trabalho
   df = DataFrame(CSV.File(file*".work"))
   dat = Matrix(df)

   n = ProgressBar(1:totlines)
   println("Plotting...")
   anim = @animate for i ∈ n
      #set_description(i, dat[i,1])
      plot(dat[i+1,2:9],dat[1,2:9], lc=:red, marker=:circle, markercolor = :red
          , label=dat[i,1],xticks=(0:5:30),yticks = yticks=(10:10:300))
      title!("Vento - "*file)
      xlabel!("velocidade [m/s]")
      ylabel!("Altura [m]")
      xlims!(0, 20)
      ylims!(10,300)
   end
   println("Creating gif image. This take a while ...")
   gif(anim, file*".gif", fps = 30)
   println("Done!")
end