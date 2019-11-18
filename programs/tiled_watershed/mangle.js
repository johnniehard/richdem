
const readline = require('readline')
const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout,
    terminal: false
})

const fc = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "properties": {nmd: {}},
            "geometry": {
                "type": "Polygon",
                "coordinates": [[]]
            }
        }
    ]
}

part = 0

rl.on('line', line => {
    switch(part){
        case 0:
            // First result: the coordinates of the polygon
            if (line === "---"){
                part += 1
                return
            }
            const [x, y] = line.trim().split(',')
            if(x === 'x') return
            fc.features[0].geometry.coordinates[0].push([+x, +y])
            break
        case 1:
            // Second result: ground cover classifications
            const [nmd_klass, count] = line.trim().split(',')
            if(nmd_klass === "nmd_klass") return
            fc.features[0].properties.nmd[nmd_klass] = +count * 4
            break
        default:
            return
    }
})


rl.on('close', () => {
    console.log(JSON.stringify(fc, null, 4))
})
