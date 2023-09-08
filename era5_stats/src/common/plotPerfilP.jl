#=
Project BRAMS-FURNAS
Author: Rodrigues, L.F. [LFR]
e-mail: luizfrodrigues@protonmail.com

This program read a CSV output of BRAMS' model, and extract the wind,
and than plot a vertical data in animated form to use in Furnas-BRAMs page

You will need the Julia Language installed, get in https://julialang.org/downloads/
and You must install the packages CSV, DAtaFrames and Plots. To install the packages
use the command using Pkg and Pkg.add("Package Name") inside Julia.

License: CC-GPL 3.0

=#
using CSV
using DataFrames
using Plots
using ArgParse
using ProgressBars
using Distributed

function justPlot(dat,i)
   plot(dat[i+1,2:9],dat[1,2:9], lc=:red, marker=:circle, markercolor = :red, label=dat[i,1]
   ,xticks=(0:5:30),yticks = yticks=(10:10:300),title = "vento", xlabel="x", ylabel = "y"
   ,xlims=(0,20),ylims=(10,300))
   return 0
end

function doOneFile(file)
   println("File input: ",file, Threads.threadid() )

   totlines = 0

   println("Adjusting header of CSV...",Threads.threadid())
   lines = readlines(file) #Le o arquivo CSV
   f = open(file*".work","w")
   iline = 0
   #Percorre as linhas do arquivo para ajustar o cabeçalho
   for line in lines
      iline = iline+1
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
      totlines = totlines+1
      # Escreve no arquivo de trabalho
      write(f,line*"\n")
   end
   close(f)
   println("Number of times: ",totlines,"...",Threads.threadid())

   # Abre o arquivo de trabalho
   df = DataFrame(CSV.File(file*".work"))
   dat = Matrix(df)

   n = (1:totlines)
   #n = ProgressBar(1:totlines)
   println("Plotting...",Threads.threadid())
   anim = @animate for i ∈ n
      println("step ",i," of ",totlines)
      x = justPlot(dat,i)
      #title!("Vento")
      #xlabel!("velocidade [m/s]")
      #ylabel!("Altura [m]")
      #xlims!(0, 20)
      #ylims!(10,300)
   end
   println("Creating gif image. This take a while ...")
   gif(anim, file*".gif", fps = 30)
   println("Done!")
end

nargs = length(ARGS)
file = ARGS

nth = Threads.nthreads()
println("Number of threads: ",nth)
println("Number of files  : ",nargs)

Threads.@threads for i = 1:nargs
   doOneFile(file[Threads.threadid()])
end
