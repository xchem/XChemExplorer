function create_view(div_name,pdb_bound,event_name,lig_name) {
// Create NGL Stage object
var stage = new NGL.Stage(div_name);
// Handle window resizing
window.addEventListener( "resize", function( event ){
    stage.handleResize();
}, false );
// Code for example: test/map-shift
Promise.all( [
    //stage.loadFile( "data://"+event_name+".ccp4" ),
    //stage.loadFile( "data://"+event_name+".pdb" )
    stage.loadFile( event_name ),
    stage.loadFile( pdb_bound ),

] ).then( function( ol ){
    var map = ol[ 0 ];
    var struc = ol[ 1 ];
    struc.autoView(lig_name)
    var surfRepr = map.addRepresentation( "surface", {
        boxSize: 10,
        useWorker: false,
        wrap: true,
        color: "skyblue",
        contour: true
    } );
    struc.addRepresentation( "cartoon" );
    struc.addRepresentation( "ball+stick", { sele: "hetero" } );
    stage.setFocus( 95 );
    stage.mouseObserver.signals.scrolled.add( function( delta ){
        if( stage.mouseObserver.altKey ){
            var d = Math.sign( delta ) / 5;
            var l = surfRepr.getParameters().isolevel;
            surfRepr.setParameters( { isolevel: l + d } );
        }
    } );
} );
};

