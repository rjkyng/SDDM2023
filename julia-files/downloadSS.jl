using CSV
using SparseArrays
using Statistics
using LinearAlgebra
using Laplacians
using MAT
using Downloads

# This function downloads the matrix `matname` from
# the suitesparse website.
# If the download is successful, it returns the matrix,
# otherwise it returns -1
# matname should be the full name of the matrix,
# see ssSddm.csv for more information.
function downloadSS(matname; remove=true)
    url = "http://sparse-files.engr.tamu.edu/mat/"
    outname = string(split(matname, '/')[end], ".mat")
    maturl = string(url, matname, ".mat")

    try
        Downloads.download(maturl, outname)
        
        file = MAT.matopen(outname)
        M = read(file, "Problem")["A"]

        if remove
            rm(outname)
        end
        
        return M
    catch
        return -1
    end
end