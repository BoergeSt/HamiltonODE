function plotOrbits(Qs,names="";kwargs...)
    X = [Q[1,:] for Q in Qs]
    Y = [Q[2,:] for Q in Qs]
    return plot(X,Y,label=hcat(names...);kwargs...)
end

function plotT(ts,Qs,titles)
    p=[]
    for (t,Q,title) in zip(ts,Qs,titles)
        push!(p,plot(t,Q',title=title))
    end
    return plot(p...)
end

plotT(ts,Qs,title::String)=plotT(ts,Qs,title,1:length(ts))
plotT(ts,Qs,title::String,extension)=plotT(ts,Qs,["$title $(extension[i])" for i in 1:length(ts)])

plotT(ts,Qs)=plotT(ts,Qs,"plot ")



println("Finished loading HamiltonPlot")