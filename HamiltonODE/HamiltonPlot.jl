function plotOrbits(Qs,names="";index=0,kwargs...)
    if index != 0
        X = [Q[1,:][index] for Q in Qs]
        Y = [Q[2,:][index] for Q in Qs]
    else
        X = [Q[1,:] for Q in Qs]
        Y = [Q[2,:] for Q in Qs]
    end
    return plot(X,Y,label=hcat(names...);kwargs...)
end

function plotT(ts,Qs,titles;kwargs...)
    p=[]
    for (t,Q,title) in zip(ts,Qs,titles)
        push!(p,plot(t,Q',title=title);kwargs...)
    end
    return plot(p...)
end

plotT(ts,Qs,title::String;kwargs...)=plotT(ts,Qs,title,1:length(ts);kwargs...)
plotT(ts,Qs,title::String,extension;kwargs...)=plotT(ts,Qs,["$title $(extension[i])" for i in 1:length(ts)];kwargs...)

plotT(ts,Qs;kwargs...)=plotT(ts,Qs,"plot ";kwargs...)



println("Finished loading HamiltonPlot")
