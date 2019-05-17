
function linearSplineInterpolation(fineT,fineX,T)
    if length(T)>1
        return map(T -> linearSplineInterpolation(fineT,fineX,T),T)
    end
    i = findfirst(fineT.>T)
    if isnothing(i)
        return fineX[end]
    end
    if i<=1
        return fineX[1]
    end
    return fineX[i-1]+(T-fineT[i-1])/(fineT[i]-fineT[i-1])*(fineX[i]-fineX[i-1])
end
