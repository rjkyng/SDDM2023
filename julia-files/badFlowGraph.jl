function speye(n)
    return sparse(I, n, n)
end


function identifyVerts(A::SparseMatrixCSC{Tv, Ti}, indLists::Array{Ti,2}) where {Tv, Ti}
    n = A.n
    P = speye(n)

    #each column of indLists is taken to be a set of verts that will be identified
    
    for k = 1:size(indLists)[2]
        inds = indLists[:,k]
        P[inds[1],inds] .= 1
    end

    vertsToRemove = indLists[2:end,:][:]
    vertsToKeep = setdiff(collect(1:n),vertsToRemove)
    P = P[vertsToKeep,:]
    return P*A*P'
end

function identifyVertsEdges(startPoints::Array{Ti,1}, endPoints::Array{Ti,1}, indLists::Array{Ti,2}) where {Tv, Ti}
    for k = 1:size(indLists)[2]
        g(v) = v in indLists[2:end][k]
        startPoints[findall(g, startPoints)] .= indLists[1][k]
    end

end

#=
function badFlowLevel(k::Int64)

    Apaths = kron(speye(k),pathGraph(k+1))
    l = k+1 #path vertex count - input to the put is n

    endIds =  [ [1;(collect(1:k-1).*l.+1)] [l;collect(2:k).*l] ]
    Akk = identifyVerts(Apaths,endIds)
    nkk = k*l-2*k+2 # k paths * len l - 2*k vertices identified with +2 vertices

    Alev = kron(speye(2),Akk)

    nlev = 2*nkk

    s1 = 1;
    t1 = l; #path ends have been identified with vertex 1 and k+1

    s2 = nkk+1;
    t2 = nkk+l;

    # add first edge between the pair of two kxk path graphs
    Alev[s1,s2] = 1;
    Alev[s2,s1] = 1;

    # add second edge
    Alev[t1,t2] = 1;
    Alev[t2,t1] = 1;

    #this is one level of the bad flow graph
    # need to connect several together

    return (Alev,s1,t1,s2,t2)
end


function badFlowGraph(k)
    A = sparse([ [0 1]; [1 0] ])
    s = 1
    t = 2

    kkEdges = [[]]
    shortEdges = [[(1,2)]]

    for i = 2:2:k
        (Alev,s1,t1,s2,t2) = badFlowLevel(i)
        nOld = size(A)[1]
        nLev = size(Alev)[1]
        Z1 = spzeros(nOld,nLev)
        Z2 = spzeros(nLev,nOld)
        #Anext = [[Alev Z1]; [Z2 A]]
        Anext = [[A Z1]; [Z2 Alev]]
        A = identifyVerts(Anext, [[s; t1+nOld] [t; s2+nOld]])
        push!(shortEdges, [(s1+nOld, t), (s, t2+nOld - 2)]) # it's -2 because we removed 2 vertices before t2+nOld
        s = s1+nOld
        t = t2+nOld - 2 # it's -2 because we removed 2 vertices before t2+nOld
    end

    return (A,s,t)
end
=#

function badFlowLevel(k::Int64)

    Apaths = kron(speye(k),pathGraph(k+1))
    l = k+1 #path vertex count - input to the put is n

    endIds =  [ [1;(collect(1:k-1).*l.+1)] [l;collect(2:k).*l] ]
    Akk = identifyVerts(Apaths,endIds)
    nkk = k*l-2*k+2 # k paths * len l - 2*k vertices identified with +2 vertices

    Alev = kron(speye(2),Akk)

    nlev = 2*nkk

    s1 = 1;
    t1 = l; #path ends have been identified with vertex 1 and k+1

    s2 = nkk+1;
    t2 = nkk+l;

    AlevKK = copy(Alev)

    # add first edge between the pair of two kxk path graphs
    Alev[s1,s2] = 1;
    Alev[s2,s1] = 1;

    # add second edge
    Alev[t1,t2] = 1;
    Alev[t2,t1] = 1;

    #this is one level of the bad flow graph
    # need to connect several together

    return (Alev, AlevKK,s1,t1,s2,t2)
end

function badFlowGraph(k::Int64)
    A = sparse([ [0 1]; [1 0] ])
    s = 1
    t = 2

    kkEdges = [[]]
    shortEdges = [[(1,2)]]

    for i = 2:2:k
        (Alev, AlevKK,s1,t1,s2,t2) = badFlowLevel(i)
        nOld = size(A)[1]
        nLev = size(Alev)[1]
        Z1 = spzeros(nOld,nLev)
        Z2 = spzeros(nLev,nOld)
        #Anext = [[Alev Z1]; [Z2 A]]
        Anext = [[A Z1]; [Z2 Alev]]
        A = identifyVerts(Anext, [[s; t1+nOld] [t; s2+nOld]])
        push!(shortEdges, [(t, s1+nOld), (s, t2+nOld - 2)]) # it's -2 because we removed 2 vertices before t2+nOld

        AnextKK = [[spzeros(nOld, nOld) Z1]; [Z2 AlevKK]]
        AKK = identifyVerts(AnextKK, [[s; t1+nOld] [t; s2+nOld]])
        rowsKK, colsKK, _ = findnz(triu(AKK))
        kkEdgesLev = []
        for x in 1:length(rowsKK)
            push!(kkEdgesLev, (rowsKK[x], colsKK[x]))
        end
        push!(kkEdges, kkEdgesLev)
        s = s1+nOld
        t = t2+nOld - 2 # it's -2 because we removed 2 vertices before t2+nOld
    end

    return (A,kkEdges,shortEdges,s,t)
end

# edges are from ai to aj
function edgeVertexMatFromEdgeList(ai::Array{Ti, 1}, aj::Array{Ti, 1}, n::Int) where {Ti, Tv}
    m = length(ai)
    return sparse(collect(1:m),aj,1.0,m,n) - sparse(collect(1:m),ai,1.0,m,n)
end


# only for spielman's graph
function fixOrientation!(A::SparseMatrixCSC{Tv, Ti}, kkEdges::Array{}, shortEdges::Array{}, s::Ti, t::Ti) where {Tv, Ti}
    n = A.n
    m = Int(length(A.nzval) / 2)
    ai = zeros(m)
    aj = zeros(m)
    maxFlow, _ = maxflow(A, s, t)
    cnt = 0
    for l in 1:length(kkEdges)
        for q in 1:length(kkEdges[l])
            cnt += 1
            (i, j) = kkEdges[l][q]
            if maxFlow[i, j] > 0
                ai[cnt] = i
                aj[cnt] = j
            else
                ai[cnt] = j
                aj[cnt] = i
                kkEdges[l][q] = (j, i)
            end
        end
    end

    for l in 1:length(shortEdges)
        for q in 1:length(shortEdges[l])
            cnt += 1
            (i, j) = shortEdges[l][q]
            if maxFlow[i, j] > 0
                ai[cnt] = i
                aj[cnt] = j
            else
                ai[cnt] = j
                aj[cnt] = i
                shortEdges[l][q] = (j, i)
            end
        end
    end

    return ai, aj
end

function R2Flow(R::Array{Tv, 1}, kkEdges::Array{}, shortEdges::Array{}, m::Ti) where {Tv, Ti}
    f = zeros(m)
    cnt = 0
    for l in 1:length(kkEdges)
        for q in 1:length(kkEdges[l])
            cnt += 1
            f[cnt] = (R[l] + R[l - 1]) / (2 * 2 * (l - 1))
        end
    end
    cnt += 1
    f[cnt] = R[1]
    for l in 2:length(shortEdges)
        for q in 1:length(shortEdges[l])
            cnt += 1
            f[cnt] = (R[l] - R[l - 1]) / 2
        end
    end
    return f
end