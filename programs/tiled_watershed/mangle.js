
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
            "properties": {},
            "geometry": {
                "type": "Polygon",
                "coordinates": [[]]
            }
        }
    ]
}

rl.on('line', line => {
    const [x, y] = line.trim().split(',')
    if(x === 'x') return
    fc.features[0].geometry.coordinates[0].push([+x, +y])
})


rl.on('close', () => {
    console.log(JSON.stringify(fc, null, 4))
})
