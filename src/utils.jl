
using HTTP
using JSON

export RESTAPI

struct RESTAPI
    baseurl
end

function (m::RESTAPI)(cmd)
    try
        url = "$(m.baseurl)/$cmd"
        response = HTTP.get(url)
        return JSON.parse(String(response.body))
    catch e
        return "error for $url occurred: $e"
    end
end
