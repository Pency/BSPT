#include "html_utils.h"

void flotr2_js(FILE* output)
{
	const char* flotr2=
"/** bean.js - copyright Jacob Thornton 2011 (MIT License) **/ /*global module:true, define:true*/\n"\

"!function(a,b,c){typeof module!=\"undefined\"?module.exports=c(a,b):typeof define==\"function\"&&typeof d"\
"efine.amd==\"object\"?define(c):b[a]=c(a,b)}(\"bean\",this,function(a,b){var c=window,d=b[a],e=/over|out/"\
",f=/[^\\.]*(?=\\..*)\\.|.*/,g=/\\..*/,h=\"addEventListener\",i=\"attachEvent\",j=\"removeEventListener\","\
"k=\"detachEvent\",l=document||{},m=l.documentElement||{},n=m[h],o=n?h:i,p=Array.prototype.slice,q=/click|"\
"mouse|menu|drag|drop/i,r=/^touch|^gesture/i,s={one:1},t=function(a,b,c){for(c=0;c<b.length;c++)a[b[c]]=1;"\
"return a}({},(\"click dblclick mouseup mousedown contextmenu mousewheel DOMMouseScroll mouseover mouseout"\
" mousemove selectstart selectend keydown keypress keyup orientationchange focus blur change reset select "\
"submit load unload beforeunload resize move DOMContentLoaded readystatechange error abort scroll \"+(n?\""\
"show input invalid touchstart touchmove touchend touchcancel gesturestart gesturechange gestureend messag"\
"e readystatechange pageshow pagehide popstate hashchange offline online afterprint beforeprint dragstart "\
"dragenter dragover dragleave drag drop dragend loadstart progress suspend emptied stalled loadmetadata lo"\
"adeddata canplay canplaythrough playing waiting seeking seeked ended durationchange timeupdate play pause"\
" ratechange volumechange cuechange checking noupdate downloading cached updateready obsolete \":\"\")).sp"\
"lit(\" \")),u=function(){function a(a,b){while((b=b.parentNode)!==null)if(b===a)return!0;return!1}functio"\
"n b(b){var c=b.relatedTarget;return c?c!==this&&c.prefix!==\"xul\"&&!/document/.test(this.toString())&&!a"\
"(this,c):c===null}return{mouseenter:{base:\"mouseover\",condition:b},mouseleave:{base:\"mouseout\",condit"\
"ion:b},mousewheel:{base:/Firefox/.test(navigator.userAgent)?\"DOMMouseScroll\":\"mousewheel\"}}}(),v=func"\
"tion(){var a=\"altKey attrChange attrName bubbles cancelable ctrlKey currentTarget detail eventPhase getM"\
"odifierState isTrusted metaKey relatedNode relatedTarget shiftKey srcElement target timeStamp type view w"\
"hich\".split(\" \"),b=a.concat(\"button buttons clientX clientY dataTransfer fromElement offsetX offsetY "\
"pageX pageY screenX screenY toElement\".split(\" \")),c=a.concat(\"char charCode key keyCode\".split(\" \"\
"")),d=a.concat(\"touches targetTouches changedTouches scale rotation\".split(\" \")),f=\"preventDefault\""\
",g=function(a){return function(){a[f]?a[f]():a.returnValue=!1}},h=\"stopPropagation\",i=function(a){retur"\
"n function(){a[h]?a[h]():a.cancelBubble=!0}},j=function(a){return function(){a[f](),a[h](),a.stopped=!0}}"\
",k=function(a,b,c){var d,e;for(d=c.length;d--;)e=c[d],!(e in b)&&e in a&&(b[e]=a[e])};return function(n,o"\
"){var p={originalEvent:n,isNative:o};if(!n)return p;var s,t=n.type,u=n.target||n.srcElement;p[f]=g(n),p[h"\
"]=i(n),p.stop=j(p),p.target=u&&u.nodeType===3?u.parentNode:u;if(o){if(t.indexOf(\"key\")!==-1)s=c,p.keyCo"\
"de=n.which||n.keyCode;else if(q.test(t)){s=b,p.rightClick=n.which===3||n.button===2,p.pos={x:0,y:0};if(n."\
"pageX||n.pageY)p.clientX=n.pageX,p.clientY=n.pageY;else if(n.clientX||n.clientY)p.clientX=n.clientX+l.bod"\
"y.scrollLeft+m.scrollLeft,p.clientY=n.clientY+l.body.scrollTop+m.scrollTop;e.test(t)&&(p.relatedTarget=n."\
"relatedTarget||n[(t===\"mouseover\"?\"from\":\"to\")+\"Element\"])}else r.test(t)&&(s=d);k(n,p,s||a)}retu"\
"rn p}}(),w=function(a,b){return!n&&!b&&(a===l||a===c)?m:a},x=function(){function a(a,b,c,d,e){this.elemen"\
"t=a,this.type=b,this.handler=c,this.original=d,this.namespaces=e,this.custom=u[b],this.isNative=t[b]&&a[o"\
"],this.eventType=n||this.isNative?b:\"propertychange\",this.customType=!n&&!this.isNative&&b,this.target="\
"w(a,this.isNative),this.eventSupport=this.target[o]}return a.prototype={inNamespaces:function(a){var b,c;"\
"if(!a)return!0;if(!this.namespaces)return!1;for(b=a.length;b--;)for(c=this.namespaces.length;c--;)if(a[b]"\
"===this.namespaces[c])return!0;return!1},matches:function(a,b,c){return this.element===a&&(!b||this.origi"\
"nal===b)&&(!c||this.handler===c)}},a}(),y=function(){var a={},b=function(c,d,e,f,g){if(!d||d===\"*\")for("\
"var h in a)h.charAt(0)===\"$\"&&b(c,h.substr(1),e,f,g);else{var i=0,j,k=a[\"$\"+d],l=c===\"*\";if(!k)retu"\
"rn;for(j=k.length;i<j;i++)if(l||k[i].matches(c,e,f))if(!g(k[i],k,i,d))return}},c=function(b,c,d){var e,f="\
"a[\"$\"+c];if(f)for(e=f.length;e--;)if(f[e].matches(b,d,null))return!0;return!1},d=function(a,c,d){var e="\
"[];return b(a,c,d,null,function(a){return e.push(a)}),e},e=function(b){return(a[\"$\"+b.type]||(a[\"$\"+b"\
".type]=[])).push(b),b},f=function(c){b(c.element,c.type,null,c.handler,function(b,c,d){return c.splice(d,"\
"1),c.length===0&&delete a[\"$\"+b.type],!1})},g=function(){var b,c=[];for(b in a)b.charAt(0)===\"$\"&&(c="\
"c.concat(a[b]));return c};return{has:c,get:d,put:e,del:f,entries:g}}(),z=n?function(a,b,c,d){a[d?h:j](b,c"\
",!1)}:function(a,b,c,d,e){e&&d&&a[\"_on\"+e]===null&&(a[\"_on\"+e]=0),a[d?i:k](\"on\"+b,c)},A=function(a,"\
"b,d){return function(e){return e=v(e||((this.ownerDocument||this.document||this).parentWindow||c).event,!"\
"0),b.apply(a,[e].concat(d))}},B=function(a,b,d,e,f,g){return function(h){if(e?e.apply(this,arguments):n?!"\
"0:h&&h.propertyName===\"_on\"+d||!h)h&&(h=v(h||((this.ownerDocument||this.document||this).parentWindow||c"\
").event,g)),b.apply(a,h&&(!f||f.length===0)?arguments:p.call(arguments,h?0:1).concat(f))}},C=function(a,b"\
",c,d,e){return function(){a(b,c,e),d.apply(this,arguments)}},D=function(a,b,c,d){var e,f,h,i=b&&b.replace"\
"(g,\"\"),j=y.get(a,i,c);for(e=0,f=j.length;e<f;e++)j[e].inNamespaces(d)&&((h=j[e]).eventSupport&&z(h.targ"\
"et,h.eventType,h.handler,!1,h.type),y.del(h))},E=function(a,b,c,d,e){var h,i=b.replace(g,\"\"),j=b.replac"\
"e(f,\"\").split(\".\");if(y.has(a,i,c))return a;i===\"unload\"&&(c=C(D,a,i,c,d)),u[i]&&(u[i].condition&&("\
"c=B(a,c,i,u[i].condition,!0)),i=u[i].base||i),h=y.put(new x(a,i,c,d,j[0]&&j)),h.handler=h.isNative?A(a,h."\
"handler,e):B(a,h.handler,i,!1,e,!1),h.eventSupport&&z(h.target,h.eventType,h.handler,!0,h.customType)},F="\
"function(a,b,c){return function(d){var e,f,g=typeof a==\"string\"?c(a,this):a;for(e=d.target;e&&e!==this;"\
"e=e.parentNode)for(f=g.length;f--;)if(g[f]===e)return b.apply(e,arguments)}},G=function(a,b,c){var d,e,h,"\
"i,j,k=D,l=b&&typeof b==\"string\";if(l&&b.indexOf(\" \")>0){b=b.split(\" \");for(j=b.length;j--;)G(a,b[j]"\
",c);return a}h=l&&b.replace(g,\"\"),h&&u[h]&&(h=u[h].type);if(!b||l){if(i=l&&b.replace(f,\"\"))i=i.split("\
"\".\");k(a,h,c,i)}else if(typeof b==\"function\")k(a,null,b);else for(d in b)b.hasOwnProperty(d)&&G(a,d,b"\
"[d]);return a},H=function(a,b,c,d,e){var f,g,h,i,j=c,k=c&&typeof c==\"string\";if(b&&!c&&typeof b==\"obje"\
"ct\")for(f in b)b.hasOwnProperty(f)&&H.apply(this,[a,f,b[f]]);else{i=arguments.length>3?p.call(arguments,"\
"3):[],g=(k?c:b).split(\" \"),k&&(c=F(b,j=d,e))&&(i=p.call(i,1)),this===s&&(c=C(G,a,b,c,j));for(h=g.length"\
";h--;)E(a,g[h],c,j,i)}return a},I=function(){return H.apply(s,arguments)},J=n?function(a,b,d){var e=l.cre"\
"ateEvent(a?\"HTMLEvents\":\"UIEvents\");e[a?\"initEvent\":\"initUIEvent\"](b,!0,!0,c,1),d.dispatchEvent(e"\
")}:function(a,b,c){c=w(c,a),a?c.fireEvent(\"on\"+b,l.createEventObject()):c[\"_on\"+b]++},K=function(a,b,"\
"c){var d,e,h,i,j,k=b.split(\" \");for(d=k.length;d--;){b=k[d].replace(g,\"\");if(i=k[d].replace(f,\"\"))i"\
"=i.split(\".\");if(!i&&!c&&a[o])J(t[b],b,a);else{j=y.get(a,b),c=[!1].concat(c);for(e=0,h=j.length;e<h;e++"\
")j[e].inNamespaces(i)&&j[e].handler.apply(a,c)}}return a},L=function(a,b,c){var d=0,e=y.get(b,c),f=e.leng"\
"th;for(;d<f;d++)e[d].original&&H(a,e[d].type,e[d].original);return a},M={add:H,one:I,remove:G,clone:L,fir"\
"e:K,noConflict:function(){return b[a]=d,this}};if(c[i]){var N=function(){var a,b=y.entries();for(a in b)b"\
"[a].type&&b[a].type!==\"unload\"&&G(b[a].element,b[a].type);c[k](\"onunload\",N),c.CollectGarbage&&c.Coll"\
"ectGarbage()};c[i](\"onunload\",N)}return M});\n"\

"/**  Underscore.js 1.1.7  (c) 2011 Jeremy Ashkenas, DocumentCloud Inc. (MIT license) **/\n"\

"(function(){var a=this,b=a._,c={},d=Array.prototype,e=Object.prototype,f=Function.prototype,g=d.slice,h=d"\
".unshift,i=e.toString,j=e.hasOwnProperty,k=d.forEach,l=d.map,m=d.reduce,n=d.reduceRight,o=d.filter,p=d.ev"\
"ery,q=d.some,r=d.indexOf,s=d.lastIndexOf,t=Array.isArray,u=Object.keys,v=f.bind,w=function(a){return new "\
"B(a)};typeof module!=\"undefined\"&&module.exports?(module.exports=w,w._=w):a._=w,w.VERSION=\"1.1.7\";var"\
" x=w.each=w.forEach=function(a,b,d){if(a==null)return;if(k&&a.forEach===k)a.forEach(b,d);else if(a.length"\
"===+a.length){for(var e=0,f=a.length;e<f;e++)if(e in a&&b.call(d,a[e],e,a)===c)return}else for(var g in a"\
")if(j.call(a,g)&&b.call(d,a[g],g,a)===c)return};w.map=function(a,b,c){var d=[];return a==null?d:l&&a.map="\
"==l?a.map(b,c):(x(a,function(a,e,f){d[d.length]=b.call(c,a,e,f)}),d)},w.reduce=w.foldl=w.inject=function("\
"a,b,c,d){var e=c!==void 0;a==null&&(a=[]);if(m&&a.reduce===m)return d&&(b=w.bind(b,d)),e?a.reduce(b,c):a."\
"reduce(b);x(a,function(a,f,g){e?c=b.call(d,c,a,f,g):(c=a,e=!0)});if(!e)throw new TypeError(\"Reduce of em"\
"pty array with no initial value\");return c},w.reduceRight=w.foldr=function(a,b,c,d){a==null&&(a=[]);if(n"\
"&&a.reduceRight===n)return d&&(b=w.bind(b,d)),c!==void 0?a.reduceRight(b,c):a.reduceRight(b);var e=(w.isA"\
"rray(a)?a.slice():w.toArray(a)).reverse();return w.reduce(e,b,c,d)},w.find=w.detect=function(a,b,c){var d"\
";return y(a,function(a,e,f){if(b.call(c,a,e,f))return d=a,!0}),d},w.filter=w.select=function(a,b,c){var d"\
"=[];return a==null?d:o&&a.filter===o?a.filter(b,c):(x(a,function(a,e,f){b.call(c,a,e,f)&&(d[d.length]=a)}"\
"),d)},w.reject=function(a,b,c){var d=[];return a==null?d:(x(a,function(a,e,f){b.call(c,a,e,f)||(d[d.lengt"\
"h]=a)}),d)},w.every=w.all=function(a,b,d){var e=!0;return a==null?e:p&&a.every===p?a.every(b,d):(x(a,func"\
"tion(a,f,g){if(!(e=e&&b.call(d,a,f,g)))return c}),e)};var y=w.some=w.any=function(a,b,d){b=b||w.identity;"\
"var e=!1;return a==null?e:q&&a.some===q?a.some(b,d):(x(a,function(a,f,g){if(e|=b.call(d,a,f,g))return c})"\
",!!e)};w.include=w.contains=function(a,b){var c=!1;return a==null?c:r&&a.indexOf===r?a.indexOf(b)!=-1:(y("\
"a,function(a){if(c=a===b)return!0}),c)},w.invoke=function(a,b){var c=g.call(arguments,2);return w.map(a,f"\
"unction(a){return(b.call?b||a:a[b]).apply(a,c)})},w.pluck=function(a,b){return w.map(a,function(a){return"\
" a[b]})},w.max=function(a,b,c){if(!b&&w.isArray(a))return Math.max.apply(Math,a);var d={computed:-Infinit"\
"y};return x(a,function(a,e,f){var g=b?b.call(c,a,e,f):a;g>=d.computed&&(d={value:a,computed:g})}),d.value"\
"},w.min=function(a,b,c){if(!b&&w.isArray(a))return Math.min.apply(Math,a);var d={computed:Infinity};retur"\
"n x(a,function(a,e,f){var g=b?b.call(c,a,e,f):a;g<d.computed&&(d={value:a,computed:g})}),d.value},w.sortB"\
"y=function(a,b,c){return w.pluck(w.map(a,function(a,d,e){return{value:a,criteria:b.call(c,a,d,e)}}).sort("\
"function(a,b){var c=a.criteria,d=b.criteria;return c<d?-1:c>d?1:0}),\"value\")},w.groupBy=function(a,b){v"\
"ar c={};return x(a,function(a,d){var e=b(a,d);(c[e]||(c[e]=[])).push(a)}),c},w.sortedIndex=function(a,b,c"\
"){c||(c=w.identity);var d=0,e=a.length;while(d<e){var f=d+e>>1;c(a[f])<c(b)?d=f+1:e=f}return d},w.toArray"\
"=function(a){return a?a.toArray?a.toArray():w.isArray(a)?g.call(a):w.isArguments(a)?g.call(a):w.values(a)"\
":[]},w.size=function(a){return w.toArray(a).length},w.first=w.head=function(a,b,c){return b!=null&&!c?g.c"\
"all(a,0,b):a[0]},w.rest=w.tail=function(a,b,c){return g.call(a,b==null||c?1:b)},w.last=function(a){return"\
" a[a.length-1]},w.compact=function(a){return w.filter(a,function(a){return!!a})},w.flatten=function(a){re"\
"turn w.reduce(a,function(a,b){return w.isArray(b)?a.concat(w.flatten(b)):(a[a.length]=b,a)},[])},w.withou"\
"t=function(a){return w.difference(a,g.call(arguments,1))},w.uniq=w.unique=function(a,b){return w.reduce(a"\
",function(a,c,d){if(0==d||(b===!0?w.last(a)!=c:!w.include(a,c)))a[a.length]=c;return a},[])},w.union=func"\
"tion(){return w.uniq(w.flatten(arguments))},w.intersection=w.intersect=function(a){var b=g.call(arguments"\
",1);return w.filter(w.uniq(a),function(a){return w.every(b,function(b){return w.indexOf(b,a)>=0})})},w.di"\
"fference=function(a,b){return w.filter(a,function(a){return!w.include(b,a)})},w.zip=function(){var a=g.ca"\
"ll(arguments),b=w.max(w.pluck(a,\"length\")),c=new Array(b);for(var d=0;d<b;d++)c[d]=w.pluck(a,\"\"+d);re"\
"turn c},w.indexOf=function(a,b,c){if(a==null)return-1;var d,e;if(c)return d=w.sortedIndex(a,b),a[d]===b?d"\
":-1;if(r&&a.indexOf===r)return a.indexOf(b);for(d=0,e=a.length;d<e;d++)if(a[d]===b)return d;return-1},w.l"\
"astIndexOf=function(a,b){if(a==null)return-1;if(s&&a.lastIndexOf===s)return a.lastIndexOf(b);var c=a.leng"\
"th;while(c--)if(a[c]===b)return c;return-1},w.range=function(a,b,c){arguments.length<=1&&(b=a||0,a=0),c=a"\
"rguments[2]||1;var d=Math.max(Math.ceil((b-a)/c),0),e=0,f=new Array(d);while(e<d)f[e++]=a,a+=c;return f},"\
"w.bind=function(a,b){if(a.bind===v&&v)return v.apply(a,g.call(arguments,1));var c=g.call(arguments,2);ret"\
"urn function(){return a.apply(b,c.concat(g.call(arguments)))}},w.bindAll=function(a){var b=g.call(argumen"\
"ts,1);return b.length==0&&(b=w.functions(a)),x(b,function(b){a[b]=w.bind(a[b],a)}),a},w.memoize=function("\
"a,b){var c={};return b||(b=w.identity),function(){var d=b.apply(this,arguments);return j.call(c,d)?c[d]:c"\
"[d]=a.apply(this,arguments)}},w.delay=function(a,b){var c=g.call(arguments,2);return setTimeout(function("\
"){return a.apply(a,c)},b)},w.defer=function(a){return w.delay.apply(w,[a,1].concat(g.call(arguments,1)))}"\
";var z=function(a,b,c){var d;return function(){var e=this,f=arguments,g=function(){d=null,a.apply(e,f)};c"\
"&&clearTimeout(d);if(c||!d)d=setTimeout(g,b)}};w.throttle=function(a,b){return z(a,b,!1)},w.debounce=func"\
"tion(a,b){return z(a,b,!0)},w.once=function(a){var b=!1,c;return function(){return b?c:(b=!0,c=a.apply(th"\
"is,arguments))}},w.wrap=function(a,b){return function(){var c=[a].concat(g.call(arguments));return b.appl"\
"y(this,c)}},w.compose=function(){var a=g.call(arguments);return function(){var b=g.call(arguments);for(va"\
"r c=a.length-1;c>=0;c--)b=[a[c].apply(this,b)];return b[0]}},w.after=function(a,b){return function(){if(-"\
"-a<1)return b.apply(this,arguments)}},w.keys=u||function(a){if(a!==Object(a))throw new TypeError(\"Invali"\
"d object\");var b=[];for(var c in a)j.call(a,c)&&(b[b.length]=c);return b},w.values=function(a){return w."\
"map(a,w.identity)},w.functions=w.methods=function(a){var b=[];for(var c in a)w.isFunction(a[c])&&b.push(c"\
");return b.sort()},w.extend=function(a){return x(g.call(arguments,1),function(b){for(var c in b)b[c]!==vo"\
"id 0&&(a[c]=b[c])}),a},w.defaults=function(a){return x(g.call(arguments,1),function(b){for(var c in b)a[c"\
"]==null&&(a[c]=b[c])}),a},w.clone=function(a){return w.isArray(a)?a.slice():w.extend({},a)},w.tap=functio"\
"n(a,b){return b(a),a},w.isEqual=function(a,b){if(a===b)return!0;var c=typeof a,d=typeof b;if(c!=d)return!"\
"1;if(a==b)return!0;if(!a&&b||a&&!b)return!1;a._chain&&(a=a._wrapped),b._chain&&(b=b._wrapped);if(a.isEqua"\
"l)return a.isEqual(b);if(b.isEqual)return b.isEqual(a);if(w.isDate(a)&&w.isDate(b))return a.getTime()===b"\
".getTime();if(w.isNaN(a)&&w.isNaN(b))return!1;if(w.isRegExp(a)&&w.isRegExp(b))return a.source===b.source&"\
"&a.global===b.global&&a.ignoreCase===b.ignoreCase&&a.multiline===b.multiline;if(c!==\"object\")return!1;i"\
"f(a.length&&a.length!==b.length)return!1;var e=w.keys(a),f=w.keys(b);if(e.length!=f.length)return!1;for(v"\
"ar g in a)if(!(g in b)||!w.isEqual(a[g],b[g]))return!1;return!0},w.isEmpty=function(a){if(w.isArray(a)||w"\
".isString(a))return a.length===0;for(var b in a)if(j.call(a,b))return!1;return!0},w.isElement=function(a)"\
"{return!!a&&a.nodeType==1},w.isArray=t||function(a){return i.call(a)===\"[object Array]\"},w.isObject=fun"\
"ction(a){return a===Object(a)},w.isArguments=function(a){return!!a&&!!j.call(a,\"callee\")},w.isFunction="\
"function(a){return!!(a&&a.constructor&&a.call&&a.apply)},w.isString=function(a){return!!(a===\"\"||a&&a.c"\
"harCodeAt&&a.substr)},w.isNumber=function(a){return!!(a===0||a&&a.toExponential&&a.toFixed)},w.isNaN=func"\
"tion(a){return a!==a},w.isBoolean=function(a){return a===!0||a===!1},w.isDate=function(a){return!!(a&&a.g"\
"etTimezoneOffset&&a.setUTCFullYear)},w.isRegExp=function(a){return!(!(a&&a.test&&a.exec)||!a.ignoreCase&&"\
"a.ignoreCase!==!1)},w.isNull=function(a){return a===null},w.isUndefined=function(a){return a===void 0},w."\
"noConflict=function(){return a._=b,this},w.identity=function(a){return a},w.times=function(a,b,c){for(var"\
" d=0;d<a;d++)b.call(c,d)},w.mixin=function(a){x(w.functions(a),function(b){D(b,w[b]=a[b])})};var A=0;w.un"\
"iqueId=function(a){var b=A++;return a?a+b:b},w.templateSettings={evaluate:/<%([\\s\\S]+?)%>/g,interpolate"\
":/<%=([\\s\\S]+?)%>/g},w.template=function(a,b){var c=w.templateSettings,d=\"var __p=[],print=function(){"\
"__p.push.apply(__p,arguments);};with(obj||{}){__p.push('\"+a.replace(/\\\\/g,\"\\\\\\\\\").replace(/'/g,\"\
""\\\\'\").replace(c.interpolate,function(a,b){return\"',\"+b.replace(/\\\\'/g,\"'\")+\",'\"}).replace(c.e"\
"valuate||null,function(a,b){return\"');\"+b.replace(/\\\\'/g,\"'\").replace(/[\\r\\n\\t]/g,\" \")+\"__p.p"\
"ush('\"}).replace(/\\r/g,\"\\\\r\").replace(/\\n/g,\"\\\\n\").replace(/\\t/g,\"\\\\t\")+\"');}return __p."\
"join('');\",e=new Function(\"obj\",d);return b?e(b):e};var B=function(a){this._wrapped=a};w.prototype=B.p"\
"rototype;var C=function(a,b){return b?w(a).chain():a},D=function(a,b){B.prototype[a]=function(){var a=g.c"\
"all(arguments);return h.call(a,this._wrapped),C(b.apply(w,a),this._chain)}};w.mixin(w),x([\"pop\",\"push\"\
"",\"reverse\",\"shift\",\"sort\",\"splice\",\"unshift\"],function(a){var b=d[a];B.prototype[a]=function()"\
"{return b.apply(this._wrapped,arguments),C(this._wrapped,this._chain)}}),x([\"concat\",\"join\",\"slice\""\
"],function(a){var b=d[a];B.prototype[a]=function(){return C(b.apply(this._wrapped,arguments),this._chain)"\
"}}),B.prototype.chain=function(){return this._chain=!0,this},B.prototype.value=function(){return this._wr"\
"apped}})();\n"\

"/** Flotr2 (c) 2012 Carl Sutherland (MIT License) **/\n"\

"(function(){var a=this,b=this.Flotr,c;c={_:_,bean:bean,isIphone:/iphone/i.test(navigator.userAgent),isIE:"\
"navigator.appVersion.indexOf(\"MSIE\")!=-1?parseFloat(navigator.appVersion.split(\"MSIE\")[1]):!1,graphTy"\
"pes:{},plugins:{},addType:function(a,b){c.graphTypes[a]=b,c.defaultOptions[a]=b.options||{},c.defaultOpti"\
"ons.defaultType=c.defaultOptions.defaultType||a},addPlugin:function(a,b){c.plugins[a]=b,c.defaultOptions["\
"a]=b.options||{}},draw:function(a,b,d,e){return e=e||c.Graph,new e(a,b,d)},merge:function(a,b){var d,e,f="\
"b||{};for(d in a)e=a[d],e&&typeof e==\"object\"?e.constructor===Array?f[d]=this._.clone(e):e.constructor!"\
"==RegExp&&!this._.isElement(e)?f[d]=c.merge(e,b?b[d]:undefined):f[d]=e:f[d]=e;return f},clone:function(a)"\
"{return c.merge(a,{})},getTickSize:function(a,b,d,e){var f=(d-b)/a,g=c.getMagnitude(f),h=10,i=f/g;return "\
"i<1.5?h=1:i<2.25?h=2:i<3?h=e===0?2:2.5:i<7.5&&(h=5),h*g},defaultTickFormatter:function(a,b){return a+\"\""\
"},defaultTrackFormatter:function(a){return\"(\"+a.x+\", \"+a.y+\")\"},engineeringNotation:function(a,b,c)"\
"{var d=[\"Y\",\"Z\",\"E\",\"P\",\"T\",\"G\",\"M\",\"k\",\"\"],e=[\"y\",\"z\",\"a\",\"f\",\"p\",\"n\",\"µ"\
"\",\"m\",\"\"],f=d.length;c=c||1e3,b=Math.pow(10,b||2);if(a===0)return 0;if(a>1)while(f--&&a>=c)a/=c;else"\
"{d=e,f=d.length;while(f--&&a<1)a*=c}return Math.round(a*b)/b+d[f]},getMagnitude:function(a){return Math.p"\
"ow(10,Math.floor(Math.log(a)/Math.LN10))},toPixel:function(a){return Math.floor(a)+.5},toRad:function(a){"\
"return-a*(Math.PI/180)},floorInBase:function(a,b){return b*Math.floor(a/b)},drawText:function(a,b,d,e,f){"\
"if(!a.fillText){a.drawText(b,d,e,f);return}f=this._.extend({size:c.defaultOptions.fontSize,color:\"#00000"\
"0\",textAlign:\"left\",textBaseline:\"bottom\",weight:1,angle:0},f),a.save(),a.translate(d,e),a.rotate(f."\
"angle),a.fillStyle=f.color,a.font=(f.weight>1?\"bold \":\"\")+f.size*1.3+\"px sans-serif\",a.textAlign=f."\
"textAlign,a.textBaseline=f.textBaseline,a.fillText(b,0,0),a.restore()},getBestTextAlign:function(a,b){ret"\
"urn b=b||{textAlign:\"center\",textBaseline:\"middle\"},a+=c.getTextAngleFromAlign(b),Math.abs(Math.cos(a"\
"))>.01&&(b.textAlign=Math.cos(a)>0?\"right\":\"left\"),Math.abs(Math.sin(a))>.01&&(b.textBaseline=Math.si"\
"n(a)>0?\"top\":\"bottom\"),b},alignTable:{\"right middle\":0,\"right top\":Math.PI/4,\"center top\":Math."\
"PI/2,\"left top\":3*(Math.PI/4),\"left middle\":Math.PI,\"left bottom\":-3*(Math.PI/4),\"center bottom\":"\
"-Math.PI/2,\"right bottom\":-Math.PI/4,\"center middle\":0},getTextAngleFromAlign:function(a){return c.al"\
"ignTable[a.textAlign+\" \"+a.textBaseline]||0},noConflict:function(){return a.Flotr=b,this}},a.Flotr=c})("\
"),Flotr.defaultOptions={colors:[\"#00A8F0\",\"#C0D800\",\"#CB4B4B\",\"#4DA74D\",\"#9440ED\"],ieBackground"\
"Color:\"#FFFFFF\",title:null,subtitle:null,shadowSize:4,defaultType:null,HtmlText:!0,fontColor:\"#545454\"\
"",fontSize:7.5,resolution:1,parseFloat:!0,xaxis:{ticks:null,minorTicks:null,showLabels:!0,showMinorLabels"\
":!1,labelsAngle:0,title:null,titleAngle:0,noTicks:5,minorTickFreq:null,tickFormatter:Flotr.defaultTickFor"\
"matter,tickDecimals:null,min:null,max:null,autoscale:!1,autoscaleMargin:0,color:null,mode:\"normal\",time"\
"Format:null,timeMode:\"UTC\",timeUnit:\"millisecond\",scaling:\"linear\",base:Math.E,titleAlign:\"center\"\
"",margin:!0},x2axis:{},yaxis:{ticks:null,minorTicks:null,showLabels:!0,showMinorLabels:!1,labelsAngle:0,t"\
"itle:null,titleAngle:90,noTicks:5,minorTickFreq:null,tickFormatter:Flotr.defaultTickFormatter,tickDecimal"\
"s:null,min:null,max:null,autoscale:!1,autoscaleMargin:0,color:null,scaling:\"linear\",base:Math.E,titleAl"\
"ign:\"center\",margin:!0},y2axis:{titleAngle:270},grid:{color:\"#545454\",backgroundColor:null,background"\
"Image:null,watermarkAlpha:.4,tickColor:\"#DDDDDD\",labelMargin:3,verticalLines:!0,minorVerticalLines:null"\
",horizontalLines:!0,minorHorizontalLines:null,outlineWidth:1,outline:\"nsew\",circular:!1},mouse:{track:!"\
"1,trackAll:!1,position:\"se\",relative:!1,trackFormatter:Flotr.defaultTrackFormatter,margin:5,lineColor:\"\
""#FF3F19\",trackDecimals:1,sensibility:2,trackY:!0,radius:3,fillColor:null,fillOpacity:.4}},function(){fu"\
"nction b(a,b,c,d){this.rgba=[\"r\",\"g\",\"b\",\"a\"];var e=4;while(-1<--e)this[this.rgba[e]]=arguments[e"\
"]||(e==3?1:0);this.normalize()}var a=Flotr._,c={aqua:[0,255,255],azure:[240,255,255],beige:[245,245,220],"\
"black:[0,0,0],blue:[0,0,255],brown:[165,42,42],cyan:[0,255,255],darkblue:[0,0,139],darkcyan:[0,139,139],d"\
"arkgrey:[169,169,169],darkgreen:[0,100,0],darkkhaki:[189,183,107],darkmagenta:[139,0,139],darkolivegreen:"\
"[85,107,47],darkorange:[255,140,0],darkorchid:[153,50,204],darkred:[139,0,0],darksalmon:[233,150,122],dar"\
"kviolet:[148,0,211],fuchsia:[255,0,255],gold:[255,215,0],green:[0,128,0],indigo:[75,0,130],khaki:[240,230"\
",140],lightblue:[173,216,230],lightcyan:[224,255,255],lightgreen:[144,238,144],lightgrey:[211,211,211],li"\
"ghtpink:[255,182,193],lightyellow:[255,255,224],lime:[0,255,0],magenta:[255,0,255],maroon:[128,0,0],navy:"\
"[0,0,128],olive:[128,128,0],orange:[255,165,0],pink:[255,192,203],purple:[128,0,128],violet:[128,0,128],r"\
"ed:[255,0,0],silver:[192,192,192],white:[255,255,255],yellow:[255,255,0]};b.prototype={scale:function(b,c"\
",d,e){var f=4;while(-1<--f)a.isUndefined(arguments[f])||(this[this.rgba[f]]*=arguments[f]);return this.no"\
"rmalize()},alpha:function(b){return!a.isUndefined(b)&&!a.isNull(b)&&(this.a=b),this.normalize()},clone:fu"\
"nction(){return new b(this.r,this.b,this.g,this.a)},limit:function(a,b,c){return Math.max(Math.min(a,c),b"\
")},normalize:function(){var a=this.limit;return this.r=a(parseInt(this.r,10),0,255),this.g=a(parseInt(thi"\
"s.g,10),0,255),this.b=a(parseInt(this.b,10),0,255),this.a=a(this.a,0,1),this},distance:function(a){if(!a)"\
"return;a=new b.parse(a);var c=0,d=3;while(-1<--d)c+=Math.abs(this[this.rgba[d]]-a[this.rgba[d]]);return c"\
"},toString:function(){return this.a>=1?\"rgb(\"+[this.r,this.g,this.b].join(\",\")+\")\":\"rgba(\"+[this."\
"r,this.g,this.b,this.a].join(\",\")+\")\"},contrast:function(){var a=1-(.299*this.r+.587*this.g+.114*this"\
".b)/255;return a<.5?\"#000000\":\"#ffffff\"}},a.extend(b,{parse:function(a){if(a instanceof b)return a;va"\
"r d;if(d=/#([a-fA-F0-9]{2})([a-fA-F0-9]{2})([a-fA-F0-9]{2})/.exec(a))return new b(parseInt(d[1],16),parse"\
"Int(d[2],16),parseInt(d[3],16));if(d=/rgb\\(\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s"\
"*\\)/.exec(a))return new b(parseInt(d[1],10),parseInt(d[2],10),parseInt(d[3],10));if(d=/#([a-fA-F0-9])([a"\
"-fA-F0-9])([a-fA-F0-9])/.exec(a))return new b(parseInt(d[1]+d[1],16),parseInt(d[2]+d[2],16),parseInt(d[3]"\
"+d[3],16));if(d=/rgba\\(\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]+(?:\\.["\
"0-9]+)?)\\s*\\)/.exec(a))return new b(parseInt(d[1],10),parseInt(d[2],10),parseInt(d[3],10),parseFloat(d["\
"4]));if(d=/rgb\\(\\s*([0-9]+(?:\\.[0-9]+)?)\\%\\s*,\\s*([0-9]+(?:\\.[0-9]+)?)\\%\\s*,\\s*([0-9]+(?:\\.[0-"\
"9]+)?)\\%\\s*\\)/.exec(a))return new b(parseFloat(d[1])*2.55,parseFloat(d[2])*2.55,parseFloat(d[3])*2.55)"\
";if(d=/rgba\\(\\s*([0-9]+(?:\\.[0-9]+)?)\\%\\s*,\\s*([0-9]+(?:\\.[0-9]+)?)\\%\\s*,\\s*([0-9]+(?:\\.[0-9]+"\
")?)\\%\\s*,\\s*([0-9]+(?:\\.[0-9]+)?)\\s*\\)/.exec(a))return new b(parseFloat(d[1])*2.55,parseFloat(d[2])"\
"*2.55,parseFloat(d[3])*2.55,parseFloat(d[4]));var e=(a+\"\").replace(/^\\s*([\\S\\s]*?)\\s*$/,\"$1\").toL"\
"owerCase();return e==\"transparent\"?new b(255,255,255,0):(d=c[e])?new b(d[0],d[1],d[2]):new b(0,0,0,0)},"\
"processColor:function(c,d){var e=d.opacity;if(!c)return\"rgba(0, 0, 0, 0)\";if(c instanceof b)return c.al"\
"pha(e).toString();if(a.isString(c))return b.parse(c).alpha(e).toString();var f=c.colors?c:{colors:c};if(!"\
"d.ctx)return a.isArray(f.colors)?b.parse(a.isArray(f.colors[0])?f.colors[0][1]:f.colors[0]).alpha(e).toSt"\
"ring():\"rgba(0, 0, 0, 0)\";f=a.extend({start:\"top\",end:\"bottom\"},f),/top/i.test(f.start)&&(d.x1=0),/"\
"left/i.test(f.start)&&(d.y1=0),/bottom/i.test(f.end)&&(d.x2=0),/right/i.test(f.end)&&(d.y2=0);var g,h,i,j"\
"=d.ctx.createLinearGradient(d.x1,d.y1,d.x2,d.y2);for(g=0;g<f.colors.length;g++)h=f.colors[g],a.isArray(h)"\
"?(i=h[0],h=h[1]):i=g/(f.colors.length-1),j.addColorStop(i,b.parse(h).alpha(e));return j}}),Flotr.Color=b}"\
"(),Flotr.Date={set:function(a,b,c,d){c=c||\"UTC\",b=\"set\"+(c===\"UTC\"?\"UTC\":\"\")+b,a[b](d)},get:fun"\
"ction(a,b,c){return c=c||\"UTC\",b=\"get\"+(c===\"UTC\"?\"UTC\":\"\")+b,a[b]()},format:function(a,b,c){fu"\
"nction f(a){return a+=\"\",a.length==1?\"0\"+a:a}if(!a)return;var d=this.get,e={h:d(a,\"Hours\",c).toStri"\
"ng(),H:f(d(a,\"Hours\",c)),M:f(d(a,\"Minutes\",c)),S:f(d(a,\"Seconds\",c)),s:d(a,\"Milliseconds\",c),d:d("\
"a,\"Date\",c).toString(),m:(d(a,\"Month\")+1).toString(),y:d(a,\"FullYear\").toString(),b:Flotr.Date.mont"\
"hNames[d(a,\"Month\",c)]},g=[],h,i=!1;for(var j=0;j<b.length;++j)h=b.charAt(j),i?(g.push(e[h]||h),i=!1):h"\
"==\"%\"?i=!0:g.push(h);return g.join(\"\")},getFormat:function(a,b){var c=Flotr.Date.timeUnits;return a<c"\
".second?\"%h:%M:%S.%s\":a<c.minute?\"%h:%M:%S\":a<c.day?b<2*c.day?\"%h:%M\":\"%b %d %h:%M\":a<c.month?\"%"\
"b %d\":a<c.year?b<c.year?\"%b\":\"%b %y\":\"%y\"},formatter:function(a,b){var c=b.options,d=Flotr.Date.ti"\
"meUnits[c.timeUnit],e=new Date(a*d);if(b.options.timeFormat)return Flotr.Date.format(e,c.timeFormat,c.tim"\
"eMode);var f=(b.max-b.min)*d,g=b.tickSize*Flotr.Date.timeUnits[b.tickUnit];return Flotr.Date.format(e,Flo"\
"tr.Date.getFormat(g,f),c.timeMode)},generator:function(a){function s(a){b(q,a,g,Flotr.floorInBase(c(q,a,g"\
"),m))}var b=this.set,c=this.get,d=this.timeUnits,e=this.spec,f=a.options,g=f.timeMode,h=d[f.timeUnit],i=a"\
".min*h,j=a.max*h,k=(j-i)/f.noTicks,l=[],m=a.tickSize,n,o,p;o=f.tickFormatter===Flotr.defaultTickFormatter"\
"?this.formatter:f.tickFormatter;for(p=0;p<e.length-1;++p){var q=e[p][0]*d[e[p][1]];if(k<(q+e[p+1][0]*d[e["\
"p+1][1]])/2&&q>=m)break}m=e[p][0],n=e[p][1],n==\"year\"&&(m=Flotr.getTickSize(f.noTicks*d.year,i,j,0),m=="\
".5&&(n=\"month\",m=6)),a.tickUnit=n,a.tickSize=m;var q=new Date(i),r=m*d[n];switch(n){case\"millisecond\""\
":s(\"Milliseconds\");break;case\"second\":s(\"Seconds\");break;case\"minute\":s(\"Minutes\");break;case\""\
"hour\":s(\"Hours\");break;case\"month\":s(\"Month\");break;case\"year\":s(\"FullYear\")}r>=d.second&&b(q,"\
"\"Milliseconds\",g,0),r>=d.minute&&b(q,\"Seconds\",g,0),r>=d.hour&&b(q,\"Minutes\",g,0),r>=d.day&&b(q,\"H"\
"ours\",g,0),r>=d.day*4&&b(q,\"Date\",g,1),r>=d.year&&b(q,\"Month\",g,0);var t=0,u=NaN,v;do{v=u,u=q.getTim"\
"e(),l.push({v:u/h,label:o(u/h,a)});if(n==\"month\")if(m<1){b(q,\"Date\",g,1);var w=q.getTime();b(q,\"Mont"\
"h\",g,c(q,\"Month\",g)+1);var x=q.getTime();q.setTime(u+t*d.hour+(x-w)*m),t=c(q,\"Hours\",g),b(q,\"Hours\"\
"",g,0)}else b(q,\"Month\",g,c(q,\"Month\",g)+m);else n==\"year\"?b(q,\"FullYear\",g,c(q,\"FullYear\",g)+m"\
"):q.setTime(u+r)}while(u<j&&u!=v);return l},timeUnits:{millisecond:1,second:1e3,minute:6e4,hour:36e5,day:"\
"864e5,month:2592e6,year:31556952e3},spec:[[1,\"millisecond\"],[20,\"millisecond\"],[50,\"millisecond\"],["\
"100,\"millisecond\"],[200,\"millisecond\"],[500,\"millisecond\"],[1,\"second\"],[2,\"second\"],[5,\"secon"\
"d\"],[10,\"second\"],[30,\"second\"],[1,\"minute\"],[2,\"minute\"],[5,\"minute\"],[10,\"minute\"],[30,\"m"\
"inute\"],[1,\"hour\"],[2,\"hour\"],[4,\"hour\"],[8,\"hour\"],[12,\"hour\"],[1,\"day\"],[2,\"day\"],[3,\"d"\
"ay\"],[.25,\"month\"],[.5,\"month\"],[1,\"month\"],[2,\"month\"],[3,\"month\"],[6,\"month\"],[1,\"year\"]"\
"],monthNames:[\"Jan\",\"Feb\",\"Mar\",\"Apr\",\"May\",\"Jun\",\"Jul\",\"Aug\",\"Sep\",\"Oct\",\"Nov\",\"D"\
"ec\"]},function(){var a=Flotr._;Flotr.DOM={addClass:function(b,c){var d=b.className?b.className:\"\";if(a"\
".include(d.split(/\\s+/g),c))return;b.className=(d?d+\" \":\"\")+c},create:function(a){return document.cr"\
"eateElement(a)},node:function(a){var b=Flotr.DOM.create(\"div\"),c;return b.innerHTML=a,c=b.children[0],b"\
".innerHTML=\"\",c},empty:function(a){a.innerHTML=\"\"},hide:function(a){Flotr.DOM.setStyles(a,{display:\""\
"none\"})},insert:function(b,c){a.isString(c)?b.innerHTML+=c:a.isElement(c)&&b.appendChild(c)},opacity:fun"\
"ction(a,b){a.style.opacity=b},position:function(a,b){return a.offsetParent?(b=this.position(a.offsetParen"\
"t),b.left+=a.offsetLeft,b.top+=a.offsetTop,b):{left:a.offsetLeft||0,top:a.offsetTop||0}},removeClass:func"\
"tion(b,c){var d=b.className?b.className:\"\";b.className=a.filter(d.split(/\\s+/g),function(a){if(a!=c)re"\
"turn!0}).join(\" \")},setStyles:function(b,c){a.each(c,function(a,c){b.style[c]=a})},show:function(a){Flo"\
"tr.DOM.setStyles(a,{display:\"\"})},size:function(a){return{height:a.offsetHeight,width:a.offsetWidth}}}}"\
"(),function(){var a=Flotr,b=a.bean;a.EventAdapter={observe:function(a,c,d){return b.add(a,c,d),this},fire"\
":function(a,c,d){return b.fire(a,c,d),typeof Prototype!=\"undefined\"&&Event.fire(a,c,d),this},stopObserv"\
"ing:function(a,c,d){return b.remove(a,c,d),this},eventPointer:function(b){if(!a._.isUndefined(b.touches)&"\
"&b.touches.length>0)return{x:b.touches[0].pageX,y:b.touches[0].pageY};if(!a._.isUndefined(b.changedTouche"\
"s)&&b.changedTouches.length>0)return{x:b.changedTouches[0].pageX,y:b.changedTouches[0].pageY};if(b.pageX|"\
"|b.pageY)return{x:b.pageX,y:b.pageY};if(b.clientX||b.clientY){var c=document,d=c.body,e=c.documentElement"\
";return{x:b.clientX+d.scrollLeft+e.scrollLeft,y:b.clientY+d.scrollTop+e.scrollTop}}}}}(),function(){var a"\
"=Flotr,b=a.DOM,c=a._,d=function(a){this.o=a};d.prototype={dimensions:function(a,b,c,d){return a?this.o.ht"\
"ml?this.html(a,this.o.element,c,d):this.canvas(a,b):{width:0,height:0}},canvas:function(b,c){if(!this.o.t"\
"extEnabled)return;c=c||{};var d=this.measureText(b,c),e=d.width,f=c.size||a.defaultOptions.fontSize,g=c.a"\
"ngle||0,h=Math.cos(g),i=Math.sin(g),j=2,k=6,l;return l={width:Math.abs(h*e)+Math.abs(i*f)+j,height:Math.a"\
"bs(i*e)+Math.abs(h*f)+k},l},html:function(a,c,d,e){var f=b.create(\"div\");return b.setStyles(f,{position"\
":\"absolute\",top:\"-10000px\"}),b.insert(f,'<div style=\"'+d+'\" class=\"'+e+' flotr-dummy-div\">'+a+\"<"\
"/div>\"),b.insert(this.o.element,f),b.size(f)},measureText:function(b,d){var e=this.o.ctx,f;return!e.fill"\
"Text||a.isIphone&&e.measure?{width:e.measure(b,d)}:(d=c.extend({size:a.defaultOptions.fontSize,weight:1,a"\
"ngle:0},d),e.save(),e.font=(d.weight>1?\"bold \":\"\")+d.size*1.3+\"px sans-serif\",f=e.measureText(b),e."\
"restore(),f)}},Flotr.Text=d}(),function(){function e(a,c,d){return b.observe.apply(this,arguments),this._"\
"handles.push(arguments),this}var a=Flotr.DOM,b=Flotr.EventAdapter,c=Flotr._,d=Flotr;Graph=function(a,e,f)"\
"{this._setEl(a),this._initMembers(),this._initPlugins(),b.fire(this.el,\"flotr:beforeinit\",[this]),this."\
"data=e,this.series=d.Series.getSeries(e),this._initOptions(f),this._initGraphTypes(),this._initCanvas(),t"\
"his._text=new d.Text({element:this.el,ctx:this.ctx,html:this.options.HtmlText,textEnabled:this.textEnable"\
"d}),b.fire(this.el,\"flotr:afterconstruct\",[this]),this._initEvents(),this.findDataRanges(),this.calcula"\
"teSpacing(),this.draw(c.bind(function(){b.fire(this.el,\"flotr:afterinit\",[this])},this))},Graph.prototy"\
"pe={destroy:function(){b.fire(this.el,\"flotr:destroy\"),c.each(this._handles,function(a){b.stopObserving"\
".apply(this,a)}),this._handles=[],this.el.graph=null},observe:e,_observe:e,processColor:function(a,b){var"\
" e={x1:0,y1:0,x2:this.plotWidth,y2:this.plotHeight,opacity:1,ctx:this.ctx};return c.extend(e,b),d.Color.p"\
"rocessColor(a,e)},findDataRanges:function(){var a=this.axes,b,e,f;c.each(this.series,function(a){f=a.getR"\
"ange(),f&&(b=a.xaxis,e=a.yaxis,b.datamin=Math.min(f.xmin,b.datamin),b.datamax=Math.max(f.xmax,b.datamax),"\
"e.datamin=Math.min(f.ymin,e.datamin),e.datamax=Math.max(f.ymax,e.datamax),b.used=b.used||f.xused,e.used=e"\
".used||f.yused)},this),!a.x.used&&!a.x2.used&&(a.x.used=!0),!a.y.used&&!a.y2.used&&(a.y.used=!0),c.each(a"\
",function(a){a.calculateRange()});var g=c.keys(d.graphTypes),h=!1;c.each(this.series,function(a){if(a.hid"\
"e)return;c.each(g,function(b){a[b]&&a[b].show&&(this.extendRange(b,a),h=!0)},this),h||this.extendRange(th"\
"is.options.defaultType,a)},this)},extendRange:function(a,b){this[a].extendRange&&this[a].extendRange(b,b."\
"data,b[a],this[a]),this[a].extendYRange&&this[a].extendYRange(b.yaxis,b.data,b[a],this[a]),this[a].extend"\
"XRange&&this[a].extendXRange(b.xaxis,b.data,b[a],this[a])},calculateSpacing:function(){var a=this.axes,b="\
"this.options,d=this.series,e=b.grid.labelMargin,f=this._text,g=a.x,h=a.x2,i=a.y,j=a.y2,k=b.grid.outlineWi"\
"dth,l,m,n,o;c.each(a,function(a){a.calculateTicks(),a.calculateTextDimensions(f,b)}),o=f.dimensions(b.tit"\
"le,{size:b.fontSize*1.5},\"font-size:1em;font-weight:bold;\",\"flotr-title\"),this.titleHeight=o.height,o"\
"=f.dimensions(b.subtitle,{size:b.fontSize},\"font-size:smaller;\",\"flotr-subtitle\"),this.subtitleHeight"\
"=o.height;for(m=0;m<b.length;++m)d[m].points.show&&(k=Math.max(k,d[m].points.radius+d[m].points.lineWidth"\
"/2));var p=this.plotOffset;g.options.margin===!1?(p.bottom=0,p.top=0):(p.bottom+=(b.grid.circular?0:g.use"\
"d&&g.options.showLabels?g.maxLabel.height+e:0)+(g.used&&g.options.title?g.titleSize.height+e:0)+k,p.top+="\
"(b.grid.circular?0:h.used&&h.options.showLabels?h.maxLabel.height+e:0)+(h.used&&h.options.title?h.titleSi"\
"ze.height+e:0)+this.subtitleHeight+this.titleHeight+k),i.options.margin===!1?(p.left=0,p.right=0):(p.left"\
"+=(b.grid.circular?0:i.used&&i.options.showLabels?i.maxLabel.width+e:0)+(i.used&&i.options.title?i.titleS"\
"ize.width+e:0)+k,p.right+=(b.grid.circular?0:j.used&&j.options.showLabels?j.maxLabel.width+e:0)+(j.used&&"\
"j.options.title?j.titleSize.width+e:0)+k),p.top=Math.floor(p.top),this.plotWidth=this.canvasWidth-p.left-"\
"p.right,this.plotHeight=this.canvasHeight-p.bottom-p.top,g.length=h.length=this.plotWidth,i.length=j.leng"\
"th=this.plotHeight,i.offset=j.offset=this.plotHeight,g.setScale(),h.setScale(),i.setScale(),j.setScale()}"\
",draw:function(a){var c=this.ctx,d;b.fire(this.el,\"flotr:beforedraw\",[this.series,this]);if(this.series"\
".length){c.save(),c.translate(this.plotOffset.left,this.plotOffset.top);for(d=0;d<this.series.length;d++)"\
"this.series[d].hide||this.drawSeries(this.series[d]);c.restore(),this.clip()}b.fire(this.el,\"flotr:after"\
"draw\",[this.series,this]),a&&a()},drawSeries:function(a){function b(a,b){var c=this.getOptions(a,b);this"\
"[b].draw(c)}var e=!1;a=a||this.series,c.each(d.graphTypes,function(c,d){a[d]&&a[d].show&&this[d]&&(e=!0,b"\
".call(this,a,d))},this),e||b.call(this,a,this.options.defaultType)},getOptions:function(a,b){var e=a[b],f"\
"=this[b],g={context:this.ctx,width:this.plotWidth,height:this.plotHeight,fontSize:this.options.fontSize,f"\
"ontColor:this.options.fontColor,textEnabled:this.textEnabled,htmlText:this.options.HtmlText,text:this._te"\
"xt,element:this.el,data:a.data,color:a.color,shadowSize:a.shadowSize,xScale:c.bind(a.xaxis.d2p,a.xaxis),y"\
"Scale:c.bind(a.yaxis.d2p,a.yaxis)};return g=d.merge(e,g),g.fillStyle=this.processColor(e.fillColor||a.col"\
"or,{opacity:e.fillOpacity}),g},getEventPosition:function(c){var d=document,e=d.body,f=d.documentElement,g"\
"=this.axes,h=this.plotOffset,i=this.lastMousePos,j=b.eventPointer(c),k=j.x-i.pageX,l=j.y-i.pageY,m,n,o;re"\
"turn\"ontouchstart\"in this.el?(m=a.position(this.overlay),n=j.x-m.left-h.left,o=j.y-m.top-h.top):(m=this"\
".overlay.getBoundingClientRect(),n=c.clientX-m.left-h.left-e.scrollLeft-f.scrollLeft,o=c.clientY-m.top-h."\
"top-e.scrollTop-f.scrollTop),{x:g.x.p2d(n),x2:g.x2.p2d(n),y:g.y.p2d(o),y2:g.y2.p2d(o),relX:n,relY:o,dX:k,"\
"dY:l,absX:j.x,absY:j.y,pageX:j.x,pageY:j.y}},clickHandler:function(a){if(this.ignoreClick)return this.ign"\
"oreClick=!1,this.ignoreClick;b.fire(this.el,\"flotr:click\",[this.getEventPosition(a),this])},mouseMoveHa"\
"ndler:function(a){if(this.mouseDownMoveHandler)return;var c=this.getEventPosition(a);b.fire(this.el,\"flo"\
"tr:mousemove\",[a,c,this]),this.lastMousePos=c},mouseDownHandler:function(a){if(this.mouseUpHandler)retur"\
"n;this.mouseUpHandler=c.bind(function(a){b.stopObserving(document,\"mouseup\",this.mouseUpHandler),b.stop"\
"Observing(document,\"mousemove\",this.mouseDownMoveHandler),this.mouseDownMoveHandler=null,this.mouseUpHa"\
"ndler=null,b.fire(this.el,\"flotr:mouseup\",[a,this])},this),this.mouseDownMoveHandler=c.bind(function(c)"\
"{var d=this.getEventPosition(c);b.fire(this.el,\"flotr:mousemove\",[a,d,this]),this.lastMousePos=d},this)"\
",b.observe(document,\"mouseup\",this.mouseUpHandler),b.observe(document,\"mousemove\",this.mouseDownMoveH"\
"andler),b.fire(this.el,\"flotr:mousedown\",[a,this]),this.ignoreClick=!1},drawTooltip:function(b,c,d,e){v"\
"ar f=this.getMouseTrack(),g=\"opacity:0.7;background-color:#000;color:#fff;display:none;position:absolute"\
";padding:2px 8px;-moz-border-radius:4px;border-radius:4px;white-space:nowrap;\",h=e.position,i=e.margin,j"\
"=this.plotOffset;c!==null&&d!==null?(e.relative?(h.charAt(0)==\"n\"?g+=\"bottom:\"+(i-j.top-d+this.canvas"\
"Height)+\"px;top:auto;\":h.charAt(0)==\"s\"&&(g+=\"top:\"+(i+j.top+d)+\"px;bottom:auto;\"),h.charAt(1)==\"\
""e\"?g+=\"left:\"+(i+j.left+c)+\"px;right:auto;\":h.charAt(1)==\"w\"&&(g+=\"right:\"+(i-j.left-c+this.can"\
"vasWidth)+\"px;left:auto;\")):(h.charAt(0)==\"n\"?g+=\"top:\"+(i+j.top)+\"px;bottom:auto;\":h.charAt(0)=="\
"\"s\"&&(g+=\"bottom:\"+(i+j.bottom)+\"px;top:auto;\"),h.charAt(1)==\"e\"?g+=\"right:\"+(i+j.right)+\"px;l"\
"eft:auto;\":h.charAt(1)==\"w\"&&(g+=\"left:\"+(i+j.left)+\"px;right:auto;\")),f.style.cssText=g,a.empty(f"\
"),a.insert(f,b),a.show(f)):a.hide(f)},clip:function(a){var b=this.plotOffset,c=this.canvasWidth,e=this.ca"\
"nvasHeight;a=a||this.ctx,d.isIE&&d.isIE<9?(a.save(),a.fillStyle=this.processColor(this.options.ieBackgrou"\
"ndColor),a.fillRect(0,0,c,b.top),a.fillRect(0,0,b.left,e),a.fillRect(0,e-b.bottom,c,b.bottom),a.fillRect("\
"c-b.right,0,b.right,e),a.restore()):(a.clearRect(0,0,c,b.top),a.clearRect(0,0,b.left,e),a.clearRect(0,e-b"\
".bottom,c,b.bottom),a.clearRect(c-b.right,0,b.right,e))},_initMembers:function(){this._handles=[],this.la"\
"stMousePos={pageX:null,pageY:null},this.plotOffset={left:0,right:0,top:0,bottom:0},this.ignoreClick=!0,th"\
"is.prevHit=null},_initGraphTypes:function(){c.each(d.graphTypes,function(a,b){this[b]=d.clone(a)},this)},"\
"_initEvents:function(){var a=this.el,d,e,f;\"ontouchstart\"in a?(d=c.bind(function(c){f=!0,b.stopObservin"\
"g(document,\"touchend\",d),b.fire(a,\"flotr:mouseup\",[event,this]),this.multitouches=null,e||this.clickH"\
"andler(c)},this),this.observe(this.overlay,\"touchstart\",c.bind(function(c){e=!1,f=!1,this.ignoreClick=!"\
"1,c.touches&&c.touches.length>1&&(this.multitouches=c.touches),b.fire(a,\"flotr:mousedown\",[event,this])"\
",this.observe(document,\"touchend\",d)},this)),this.observe(this.overlay,\"touchmove\",c.bind(function(c)"\
"{var d=this.getEventPosition(c);c.preventDefault(),e=!0,this.multitouches||c.touches&&c.touches.length>1?"\
"this.multitouches=c.touches:f||b.fire(a,\"flotr:mousemove\",[event,d,this]),this.lastMousePos=d},this))):"\
"this.observe(this.overlay,\"mousedown\",c.bind(this.mouseDownHandler,this)).observe(a,\"mousemove\",c.bin"\
"d(this.mouseMoveHandler,this)).observe(this.overlay,\"click\",c.bind(this.clickHandler,this)).observe(a,\"\
""mouseout\",function(){b.fire(a,\"flotr:mouseout\")})},_initCanvas:function(){function k(e,f){return e||("\
"e=a.create(\"canvas\"),typeof FlashCanvas!=\"undefined\"&&typeof e.getContext==\"function\"&&FlashCanvas."\
"initElement(e),e.className=\"flotr-\"+f,e.style.cssText=\"position:absolute;left:0px;top:0px;\",a.insert("\
"b,e)),c.each(i,function(b,c){a.show(e);if(f==\"canvas\"&&e.getAttribute(c)===b)return;e.setAttribute(c,b*"\
"d.resolution),e.style[c]=b+\"px\"}),e.context_=null,e}function l(a){window.G_vmlCanvasManager&&window.G_v"\
"mlCanvasManager.initElement(a);var b=a.getContext(\"2d\");return window.G_vmlCanvasManager||b.scale(d.res"\
"olution,d.resolution),b}var b=this.el,d=this.options,e=b.children,f=[],g,h,i,j;for(h=e.length;h--;)g=e[h]"\
",!this.canvas&&g.className===\"flotr-canvas\"?this.canvas=g:!this.overlay&&g.className===\"flotr-overlay\"\
""?this.overlay=g:f.push(g);for(h=f.length;h--;)b.removeChild(f[h]);a.setStyles(b,{position:\"relative\"})"\
",i={},i.width=b.clientWidth,i.height=b.clientHeight;if(i.width<=0||i.height<=0||d.resolution<=0)throw\"In"\
"valid dimensions for plot, width = \"+i.width+\", height = \"+i.height+\", resolution = \"+d.resolution;t"\
"his.canvas=k(this.canvas,\"canvas\"),this.overlay=k(this.overlay,\"overlay\"),this.ctx=l(this.canvas),thi"\
"s.ctx.clearRect(0,0,this.canvas.width,this.canvas.height),this.octx=l(this.overlay),this.octx.clearRect(0"\
",0,this.overlay.width,this.overlay.height),this.canvasHeight=i.height,this.canvasWidth=i.width,this.textE"\
"nabled=!!this.ctx.drawText||!!this.ctx.fillText},_initPlugins:function(){c.each(d.plugins,function(a,b){c"\
".each(a.callbacks,function(a,b){this.observe(this.el,b,c.bind(a,this))},this),this[b]=d.clone(a),c.each(t"\
"his[b],function(a,d){c.isFunction(a)&&(this[b][d]=c.bind(a,this))},this)},this)},_initOptions:function(a)"\
"{var e=d.clone(d.defaultOptions);e.x2axis=c.extend(c.clone(e.xaxis),e.x2axis),e.y2axis=c.extend(c.clone(e"\
".yaxis),e.y2axis),this.options=d.merge(a||{},e),this.options.grid.minorVerticalLines===null&&this.options"\
".xaxis.scaling===\"logarithmic\"&&(this.options.grid.minorVerticalLines=!0),this.options.grid.minorHorizo"\
"ntalLines===null&&this.options.yaxis.scaling===\"logarithmic\"&&(this.options.grid.minorHorizontalLines=!"\
"0),b.fire(this.el,\"flotr:afterinitoptions\",[this]),this.axes=d.Axis.getAxes(this.options);var f=[],g=[]"\
",h=this.series.length,i=this.series.length,j=this.options.colors,k=[],l=0,m,n,o,p;for(n=i-1;n>-1;--n)m=th"\
"is.series[n].color,m&&(--i,c.isNumber(m)?f.push(m):k.push(d.Color.parse(m)));for(n=f.length-1;n>-1;--n)i="\
"Math.max(i,f[n]+1);for(n=0;g.length<i;){m=j.length==n?new d.Color(100,100,100):d.Color.parse(j[n]);var q="\
"l%2==1?-1:1,r=1+q*Math.ceil(l/2)*.2;m.scale(r,r,r),g.push(m),++n>=j.length&&(n=0,++l)}for(n=0,o=0;n<h;++n"\
"){p=this.series[n],p.color?c.isNumber(p.color)&&(p.color=g[p.color].toString()):p.color=g[o++].toString()"\
",p.xaxis||(p.xaxis=this.axes.x),p.xaxis==1?p.xaxis=this.axes.x:p.xaxis==2&&(p.xaxis=this.axes.x2),p.yaxis"\
"||(p.yaxis=this.axes.y),p.yaxis==1?p.yaxis=this.axes.y:p.yaxis==2&&(p.yaxis=this.axes.y2);for(var s in d."\
"graphTypes)p[s]=c.extend(c.clone(this.options[s]),p[s]);p.mouse=c.extend(c.clone(this.options.mouse),p.mo"\
"use),c.isUndefined(p.shadowSize)&&(p.shadowSize=this.options.shadowSize)}},_setEl:function(a){if(!a)throw"\
"\"The target container doesn't exist\";if(a.graph instanceof Graph)a.graph.destroy();else if(!a.clientWid"\
"th)throw\"The target container must be visible\";a.graph=this,this.el=a}},Flotr.Graph=Graph}(),function()"\
"{function c(b){this.orientation=1,this.offset=0,this.datamin=Number.MAX_VALUE,this.datamax=-Number.MAX_VA"\
"LUE,a.extend(this,b),this._setTranslations()}function d(a){return this.offset+this.orientation*(a-this.mi"\
"n)*this.scale}function e(a){return(this.offset+this.orientation*a)/this.scale+this.min}function f(a){retu"\
"rn this.offset+this.orientation*(h(a,this.options.base)-h(this.min,this.options.base))*this.scale}functio"\
"n g(a){return j((this.offset+this.orientation*a)/this.scale+h(this.min,this.options.base),this.options.ba"\
"se)}function h(a,b){return a=Math.log(Math.max(a,Number.MIN_VALUE)),b!==Math.E&&(a/=Math.log(b)),a}functi"\
"on j(a,b){return b===Math.E?Math.exp(a):Math.pow(b,a)}var a=Flotr._,b=\"logarithmic\";c.prototype={setSca"\
"le:function(){var a=this.length;this.options.scaling==b?this.scale=a/(h(this.max,this.options.base)-h(thi"\
"s.min,this.options.base)):this.scale=a/(this.max-this.min)},calculateTicks:function(){var b=this.options;"\
"this.ticks=[],this.minorTicks=[],b.ticks?(this._cleanUserTicks(b.ticks,this.ticks),this._cleanUserTicks(b"\
".minorTicks||[],this.minorTicks)):b.mode==\"time\"?this._calculateTimeTicks():b.scaling===\"logarithmic\""\
"?this._calculateLogTicks():this._calculateTicks(),a.each(this.ticks,function(a){a.label+=\"\"}),a.each(th"\
"is.minorTicks,function(a){a.label+=\"\"})},calculateRange:function(){if(!this.used)return;var a=this,b=a."\
"options,c=b.min!==null?b.min:a.datamin,d=b.max!==null?b.max:a.datamax,e=b.autoscaleMargin;b.scaling==\"lo"\
"garithmic\"&&(c<=0&&(c=a.datamin),d<=0&&(d=c));if(d==c){var f=d?.01:1;b.min===null&&(c-=f),b.max===null&&"\
"(d+=f)}if(b.scaling===\"logarithmic\"){c<0&&(c=d/b.base);var g=Math.log(d);b.base!=Math.E&&(g/=Math.log(b"\
".base)),g=Math.ceil(g);var h=Math.log(c);b.base!=Math.E&&(h/=Math.log(b.base)),h=Math.ceil(h),a.tickSize="\
"Flotr.getTickSize(b.noTicks,h,g,b.tickDecimals===null?0:b.tickDecimals),b.minorTickFreq===null&&(g-h>10?b"\
".minorTickFreq=0:g-h>5?b.minorTickFreq=2:b.minorTickFreq=5)}else a.tickSize=Flotr.getTickSize(b.noTicks,c"\
",d,b.tickDecimals);a.min=c,a.max=d,b.min===null&&b.autoscale&&(a.min-=a.tickSize*e,a.min<0&&a.datamin>=0&"\
"&(a.min=0),a.min=a.tickSize*Math.floor(a.min/a.tickSize)),b.max===null&&b.autoscale&&(a.max+=a.tickSize*e"\
",a.max>0&&a.datamax<=0&&a.datamax!=a.datamin&&(a.max=0),a.max=a.tickSize*Math.ceil(a.max/a.tickSize)),a.m"\
"in==a.max&&(a.max=a.min+1)},calculateTextDimensions:function(a,b){var c=\"\",d,e;if(this.options.showLabe"\
"ls)for(e=0;e<this.ticks.length;++e)d=this.ticks[e].label.length,d>c.length&&(c=this.ticks[e].label);this."\
"maxLabel=a.dimensions(c,{size:b.fontSize,angle:Flotr.toRad(this.options.labelsAngle)},\"font-size:smaller"\
";\",\"flotr-grid-label\"),this.titleSize=a.dimensions(this.options.title,{size:b.fontSize*1.2,angle:Flotr"\
".toRad(this.options.titleAngle)},\"font-weight:bold;\",\"flotr-axis-title\")},_cleanUserTicks:function(b,"\
"c){var d=this,e=this.options,f,g,h,i;a.isFunction(b)&&(b=b({min:d.min,max:d.max}));for(g=0;g<b.length;++g"\
")i=b[g],typeof i==\"object\"?(f=i[0],h=i.length>1?i[1]:e.tickFormatter(f,{min:d.min,max:d.max})):(f=i,h=e"\
".tickFormatter(f,{min:this.min,max:this.max})),c[g]={v:f,label:h}},_calculateTimeTicks:function(){this.ti"\
"cks=Flotr.Date.generator(this)},_calculateLogTicks:function(){var a=this,b=a.options,c,d,e=Math.log(a.max"\
");b.base!=Math.E&&(e/=Math.log(b.base)),e=Math.ceil(e);var f=Math.log(a.min);b.base!=Math.E&&(f/=Math.log"\
"(b.base)),f=Math.ceil(f);for(i=f;i<e;i+=a.tickSize){d=b.base==Math.E?Math.exp(i):Math.pow(b.base,i);var g"\
"=d*(b.base==Math.E?Math.exp(a.tickSize):Math.pow(b.base,a.tickSize)),h=(g-d)/b.minorTickFreq;a.ticks.push"\
"({v:d,label:b.tickFormatter(d,{min:a.min,max:a.max})});for(c=d+h;c<g;c+=h)a.minorTicks.push({v:c,label:b."\
"tickFormatter(c,{min:a.min,max:a.max})})}d=b.base==Math.E?Math.exp(i):Math.pow(b.base,i),a.ticks.push({v:"\
"d,label:b.tickFormatter(d,{min:a.min,max:a.max})})},_calculateTicks:function(){var a=this,b=a.options,c=a"\
".tickSize,d=a.min,e=a.max,f=c*Math.ceil(d/c),g,h,i,j,k,l;b.minorTickFreq&&(h=c/b.minorTickFreq);for(k=0;("\
"i=j=f+k*c)<=e;++k){g=b.tickDecimals,g===null&&(g=1-Math.floor(Math.log(c)/Math.LN10)),g<0&&(g=0),i=i.toFi"\
"xed(g),a.ticks.push({v:i,label:b.tickFormatter(i,{min:a.min,max:a.max})});if(b.minorTickFreq)for(l=0;l<b."\
"minorTickFreq&&k*c+l*h<e;++l)i=j+l*h,a.minorTicks.push({v:i,label:b.tickFormatter(i,{min:a.min,max:a.max}"\
")})}},_setTranslations:function(a){this.d2p=a?f:d,this.p2d=a?g:e}},a.extend(c,{getAxes:function(a){return"\
"{x:new c({options:a.xaxis,n:1,length:this.plotWidth}),x2:new c({options:a.x2axis,n:2,length:this.plotWidt"\
"h}),y:new c({options:a.yaxis,n:1,length:this.plotHeight,offset:this.plotHeight,orientation:-1}),y2:new c("\
"{options:a.y2axis,n:2,length:this.plotHeight,offset:this.plotHeight,orientation:-1})}}}),Flotr.Axis=c}(),"\
"function(){function b(b){a.extend(this,b)}var a=Flotr._;b.prototype={getRange:function(){var a=this.data,"\
"b=a.length,c=Number.MAX_VALUE,d=Number.MAX_VALUE,e=-Number.MAX_VALUE,f=-Number.MAX_VALUE,g=!1,h=!1,i,j,k;"\
"if(b<0||this.hide)return!1;for(k=0;k<b;k++)i=a[k][0],j=a[k][1],i<c&&(c=i,g=!0),i>e&&(e=i,g=!0),j<d&&(d=j,"\
"h=!0),j>f&&(f=j,h=!0);return{xmin:c,xmax:e,ymin:d,ymax:f,xused:g,yused:h}}},a.extend(b,{getSeries:functio"\
"n(c){return a.map(c,function(c){var d;return c.data?(d=new b,a.extend(d,c)):d=new b({data:c}),d})}}),Flot"\
"r.Series=b}(),Flotr.addType(\"lines\",{options:{show:!1,lineWidth:2,fill:!1,fillBorder:!1,fillColor:null,"\
"fillOpacity:.4,steps:!1,stacked:!1},stack:{values:[]},draw:function(a){var b=a.context,c=a.lineWidth,d=a."\
"shadowSize,e;b.save(),b.lineJoin=\"round\",d&&(b.lineWidth=d/2,e=c/2+b.lineWidth/2,b.strokeStyle=\"rgba(0"\
",0,0,0.1)\",this.plot(a,e+d/2,!1),b.strokeStyle=\"rgba(0,0,0,0.2)\",this.plot(a,e,!1)),b.lineWidth=c,b.st"\
"rokeStyle=a.color,this.plot(a,0,!0),b.restore()},plot:function(a,b,c){function w(){!b&&a.fill&&o&&(p=g(o["\
"0]),d.fillStyle=a.fillStyle,d.lineTo(q,n),d.lineTo(p,n),d.lineTo(p,h(o[1])),d.fill(),a.fillBorder&&d.stro"\
"ke())}var d=a.context,e=a.width,f=a.height,g=a.xScale,h=a.yScale,i=a.data,j=a.stacked?this.stack:!1,k=i.l"\
"ength-1,l=null,m=null,n=h(0),o=null,p,q,r,s,t,u,v;if(k<1)return;d.beginPath();for(v=0;v<k;++v){if(i[v][1]"\
"===null||i[v+1][1]===null){a.fill&&v>0&&i[v][1]&&(d.stroke(),w(),o=null,d.closePath(),d.beginPath());cont"\
"inue}p=g(i[v][0]),q=g(i[v+1][0]),o===null&&(o=i[v]),j?(t=j.values[i[v][0]]||0,u=j.values[i[v+1][0]]||j.va"\
"lues[i[v][0]]||0,r=h(i[v][1]+t),s=h(i[v+1][1]+u),c&&(j.values[i[v][0]]=i[v][1]+t,v==k-1&&(j.values[i[v+1]"\
"[0]]=i[v+1][1]+u))):(r=h(i[v][1]),s=h(i[v+1][1]));if(r>f&&s>f||r<0&&s<0||p<0&&q<0||p>e&&q>e)continue;(l!="\
"p||m!=r+b)&&d.moveTo(p,r+b),l=q,m=s+b,a.steps?(d.lineTo(l+b/2,r+b),d.lineTo(l+b/2,m)):d.lineTo(l,m)}(!a.f"\
"ill||a.fill&&!a.fillBorder)&&d.stroke(),w(),d.closePath()},extendYRange:function(a,b,c,d){var e=a.options"\
";if(c.stacked&&(!e.max&&e.max!==0||!e.min&&e.min!==0)){var f=a.max,g=a.min,h=d.positiveSums||{},i=d.negat"\
"iveSums||{},j,k;for(k=0;k<b.length;k++)j=b[k][0]+\"\",b[k][1]>0?(h[j]=(h[j]||0)+b[k][1],f=Math.max(f,h[j]"\
")):(i[j]=(i[j]||0)+b[k][1],g=Math.min(g,i[j]));d.negativeSums=i,d.positiveSums=h,a.max=f,a.min=g}c.steps&"\
"&(this.hit=function(a){var b=a.data,c=a.args,d=a.yScale,e=c[0],f=b.length,g=c[1],h=e.x,i=e.relY,j;for(j=0"\
";j<f-1;j++)if(h>=b[j][0]&&h<=b[j+1][0]){Math.abs(d(b[j][1])-i)<8&&(g.x=b[j][0],g.y=b[j][1],g.index=j,g.se"\
"riesIndex=a.index);break}},this.drawHit=function(a){var b=a.context,c=a.args,d=a.data,e=a.xScale,f=c.inde"\
"x,g=e(c.x),h=a.yScale(c.y),i;d.length-1>f&&(i=a.xScale(d[f+1][0]),b.save(),b.strokeStyle=a.color,b.lineWi"\
"dth=a.lineWidth,b.beginPath(),b.moveTo(g,h),b.lineTo(i,h),b.stroke(),b.closePath(),b.restore())},this.cle"\
"arHit=function(a){var b=a.context,c=a.args,d=a.data,e=a.xScale,f=a.lineWidth,g=c.index,h=e(c.x),i=a.yScal"\
"e(c.y),j;d.length-1>g&&(j=a.xScale(d[g+1][0]),b.clearRect(h-f,i-f,j-h+2*f,2*f))})}}),Flotr.addType(\"bars"\
"\",{options:{show:!1,lineWidth:2,barWidth:1,fill:!0,fillColor:null,fillOpacity:.4,horizontal:!1,stacked:!"\
"1,centered:!0,topPadding:.1,grouped:!1},stack:{positive:[],negative:[],_positive:[],_negative:[]},draw:fu"\
"nction(a){var b=a.context;this.current+=1,b.save(),b.lineJoin=\"miter\",b.lineWidth=a.lineWidth,b.strokeS"\
"tyle=a.color,a.fill&&(b.fillStyle=a.fillStyle),this.plot(a),b.restore()},plot:function(a){var b=a.data,c="\
"a.context,d=a.shadowSize,e,f,g,h,i,j;if(b.length<1)return;this.translate(c,a.horizontal);for(e=0;e<b.leng"\
"th;e++){f=this.getBarGeometry(b[e][0],b[e][1],a);if(f===null)continue;g=f.left,h=f.top,i=f.width,j=f.heig"\
"ht,a.fill&&c.fillRect(g,h,i,j),d&&(c.save(),c.fillStyle=\"rgba(0,0,0,0.05)\",c.fillRect(g+d,h+d,i,j),c.re"\
"store()),a.lineWidth&&c.strokeRect(g,h,i,j)}},translate:function(a,b){b&&(a.rotate(-Math.PI/2),a.scale(-1"\
",1))},getBarGeometry:function(a,b,c){var d=c.horizontal,e=c.barWidth,f=c.centered,g=c.stacked?this.stack:"\
"!1,h=c.lineWidth,i=f?e/2:0,j=d?c.yScale:c.xScale,k=d?c.xScale:c.yScale,l=d?b:a,m=d?a:b,n=0,o,p,q,r,s;retu"\
"rn c.grouped&&(this.current/this.groups,l-=i,e/=this.groups,i=e/2,l=l+e*this.current-i),g&&(o=m>0?g.posit"\
"ive:g.negative,n=o[l]||n,o[l]=n+m),p=j(l-i),q=j(l+e-i),r=k(m+n),s=k(n),s<0&&(s=0),a===null||b===null?null"\
":{x:l,y:m,xScale:j,yScale:k,top:r,left:Math.min(p,q)-h/2,width:Math.abs(q-p)-h,height:s-r}},hit:function("\
"a){var b=a.data,c=a.args,d=c[0],e=c[1],f=d.x,g=d.y,h=this.getBarGeometry(f,g,a),i=h.width/2,j=h.left,k,l;"\
"for(l=b.length;l--;)k=this.getBarGeometry(b[l][0],b[l][1],a),k.y>h.y&&Math.abs(j-k.left)<i&&(e.x=b[l][0],"\
"e.y=b[l][1],e.index=l,e.seriesIndex=a.index)},drawHit:function(a){var b=a.context,c=a.args,d=this.getBarG"\
"eometry(c.x,c.y,a),e=d.left,f=d.top,g=d.width,h=d.height;b.save(),b.strokeStyle=a.color,b.lineWidth=a.lin"\
"eWidth,this.translate(b,a.horizontal),b.beginPath(),b.moveTo(e,f+h),b.lineTo(e,f),b.lineTo(e+g,f),b.lineT"\
"o(e+g,f+h),a.fill&&(b.fillStyle=a.fillStyle,b.fill()),b.stroke(),b.closePath(),b.restore()},clearHit:func"\
"tion(a){var b=a.context,c=a.args,d=this.getBarGeometry(c.x,c.y,a),e=d.left,f=d.width,g=d.top,h=d.height,i"\
"=2*a.lineWidth;b.save(),this.translate(b,a.horizontal),b.clearRect(e-i,Math.min(g,g+h)-i,f+2*i,Math.abs(h"\
")+2*i),b.restore()},extendXRange:function(a,b,c,d){this._extendRange(a,b,c,d),this.groups=this.groups+1||"\
"1,this.current=0},extendYRange:function(a,b,c,d){this._extendRange(a,b,c,d)},_extendRange:function(a,b,c,"\
"d){var e=a.options.max;if(_.isNumber(e)||_.isString(e))return;var f=a.min,g=a.max,h=c.horizontal,i=a.orie"\
"ntation,j=this.positiveSums||{},k=this.negativeSums||{},l,m,n,o;(i==1&&!h||i==-1&&h)&&c.centered&&(g=Math"\
".max(a.datamax+c.barWidth,g),f=Math.min(a.datamin-c.barWidth,f));if(c.stacked&&(i==1&&h||i==-1&&!h))for(o"\
"=b.length;o--;)l=b[o][i==1?1:0]+\"\",m=b[o][i==1?0:1],m>0?(j[l]=(j[l]||0)+m,g=Math.max(g,j[l])):(k[l]=(k["\
"l]||0)+m,f=Math.min(f,k[l]));(i==1&&h||i==-1&&!h)&&c.topPadding&&(a.max===a.datamax||c.stacked&&this.stac"\
"kMax!==g)&&(g+=c.topPadding*(g-f)),this.stackMin=f,this.stackMax=g,this.negativeSums=k,this.positiveSums="\
"j,a.max=g,a.min=f}}),Flotr.addType(\"bubbles\",{options:{show:!1,lineWidth:2,fill:!0,fillOpacity:.4,baseR"\
"adius:2},draw:function(a){var b=a.context,c=a.shadowSize;b.save(),b.lineWidth=a.lineWidth,b.fillStyle=\"r"\
"gba(0,0,0,0.05)\",b.strokeStyle=\"rgba(0,0,0,0.05)\",this.plot(a,c/2),b.strokeStyle=\"rgba(0,0,0,0.1)\",t"\
"his.plot(a,c/4),b.strokeStyle=a.color,b.fillStyle=a.fillStyle,this.plot(a),b.restore()},plot:function(a,b"\
"){var c=a.data,d=a.context,e,f,g,h,i;b=b||0;for(f=0;f<c.length;++f)e=this.getGeometry(c[f],a),d.beginPath"\
"(),d.arc(e.x+b,e.y+b,e.z,0,2*Math.PI,!0),d.stroke(),a.fill&&d.fill(),d.closePath()},getGeometry:function("\
"a,b){return{x:b.xScale(a[0]),y:b.yScale(a[1]),z:a[2]*b.baseRadius}},hit:function(a){var b=a.data,c=a.args"\
",d=c[0],e=c[1],f=d.x,g=d.y,h,j,k,l;e.best=e.best||Number.MAX_VALUE;for(i=b.length;i--;)j=this.getGeometry"\
"(b[i],a),k=j.x-a.xScale(f),l=j.y-a.yScale(g),h=Math.sqrt(k*k+l*l),h<j.z&&j.z<e.best&&(e.x=b[i][0],e.y=b[i"\
"][1],e.index=i,e.seriesIndex=a.index,e.best=j.z)},drawHit:function(a){var b=a.context,c=this.getGeometry("\
"a.data[a.args.index],a);b.save(),b.lineWidth=a.lineWidth,b.fillStyle=a.fillStyle,b.strokeStyle=a.color,b."\
"beginPath(),b.arc(c.x,c.y,c.z,0,2*Math.PI,!0),b.fill(),b.stroke(),b.closePath(),b.restore()},clearHit:fun"\
"ction(a){var b=a.context,c=this.getGeometry(a.data[a.args.index],a),d=c.z+a.lineWidth;b.save(),b.clearRec"\
"t(c.x-d,c.y-d,2*d,2*d),b.restore()}}),Flotr.addType(\"candles\",{options:{show:!1,lineWidth:1,wickLineWid"\
"th:1,candleWidth:.6,fill:!0,upFillColor:\"#00A8F0\",downFillColor:\"#CB4B4B\",fillOpacity:.5,barcharts:!1"\
"},draw:function(a){var b=a.context;b.save(),b.lineJoin=\"miter\",b.lineCap=\"butt\",b.lineWidth=a.wickLin"\
"eWidth||a.lineWidth,this.plot(a),b.restore()},plot:function(a){var b=a.data,c=a.context,d=a.xScale,e=a.yS"\
"cale,f=a.candleWidth/2,g=a.shadowSize,h=a.lineWidth,i=a.wickLineWidth,j=i%2/2,k,l,m,n,o,p,q,r,s,t,u,v,w,x"\
",y;if(b.length<1)return;for(y=0;y<b.length;y++){l=b[y],m=l[0],o=l[1],p=l[2],q=l[3],r=l[4],s=d(m-f),t=d(m+"\
"f),u=e(q),v=e(p),w=e(Math.min(o,r)),x=e(Math.max(o,r)),k=a[o>r?\"downFillColor\":\"upFillColor\"],a.fill&"\
"&!a.barcharts&&(c.fillStyle=\"rgba(0,0,0,0.05)\",c.fillRect(s+g,x+g,t-s,w-x),c.save(),c.globalAlpha=a.fil"\
"lOpacity,c.fillStyle=k,c.fillRect(s,x+h,t-s,w-x),c.restore());if(h||i)m=Math.floor((s+t)/2)+j,c.strokeSty"\
"le=k,c.beginPath(),a.barcharts?(c.moveTo(m,Math.floor(v+f)),c.lineTo(m,Math.floor(u+f)),n=Math.floor(o+f)"\
"+.5,c.moveTo(Math.floor(s)+j,n),c.lineTo(m,n),n=Math.floor(r+f)+.5,c.moveTo(Math.floor(t)+j,n),c.lineTo(m"\
",n)):(c.strokeRect(s,x+h,t-s,w-x),c.moveTo(m,Math.floor(x+h)),c.lineTo(m,Math.floor(v+h)),c.moveTo(m,Math"\
".floor(w+h)),c.lineTo(m,Math.floor(u+h))),c.closePath(),c.stroke()}},extendXRange:function(a,b,c){a.optio"\
"ns.max===null&&(a.max=Math.max(a.datamax+.5,a.max),a.min=Math.min(a.datamin-.5,a.min))}}),Flotr.addType(\"\
""gantt\",{options:{show:!1,lineWidth:2,barWidth:1,fill:!0,fillColor:null,fillOpacity:.4,centered:!0},draw"\
":function(a){var b=this.ctx,c=a.gantt.barWidth,d=Math.min(a.gantt.lineWidth,c);b.save(),b.translate(this."\
"plotOffset.left,this.plotOffset.top),b.lineJoin=\"miter\",b.lineWidth=d,b.strokeStyle=a.color,b.save(),th"\
"is.gantt.plotShadows(a,c,0,a.gantt.fill),b.restore();if(a.gantt.fill){var e=a.gantt.fillColor||a.color;b."\
"fillStyle=this.processColor(e,{opacity:a.gantt.fillOpacity})}this.gantt.plot(a,c,0,a.gantt.fill),b.restor"\
"e()},plot:function(a,b,c,d){var e=a.data;if(e.length<1)return;var f=a.xaxis,g=a.yaxis,h=this.ctx,i;for(i="\
"0;i<e.length;i++){var j=e[i][0],k=e[i][1],l=e[i][2],m=!0,n=!0,o=!0;if(k===null||l===null)continue;var p=k"\
",q=k+l,r=j-(a.gantt.centered?b/2:0),s=j+b-(a.gantt.centered?b/2:0);if(q<f.min||p>f.max||s<g.min||r>g.max)"\
"continue;p<f.min&&(p=f.min,m=!1),q>f.max&&(q=f.max,f.lastSerie!=a&&(n=!1)),r<g.min&&(r=g.min),s>g.max&&(s"\
"=g.max,g.lastSerie!=a&&(n=!1)),d&&(h.beginPath(),h.moveTo(f.d2p(p),g.d2p(r)+c),h.lineTo(f.d2p(p),g.d2p(s)"\
"+c),h.lineTo(f.d2p(q),g.d2p(s)+c),h.lineTo(f.d2p(q),g.d2p(r)+c),h.fill(),h.closePath()),a.gantt.lineWidth"\
"&&(m||o||n)&&(h.beginPath(),h.moveTo(f.d2p(p),g.d2p(r)+c),h[m?\"lineTo\":\"moveTo\"](f.d2p(p),g.d2p(s)+c)"\
",h[n?\"lineTo\":\"moveTo\"](f.d2p(q),g.d2p(s)+c),h[o?\"lineTo\":\"moveTo\"](f.d2p(q),g.d2p(r)+c),h.stroke"\
"(),h.closePath())}},plotShadows:function(a,b,c){var d=a.data;if(d.length<1)return;var e,f,g,h,i=a.xaxis,j"\
"=a.yaxis,k=this.ctx,l=this.options.shadowSize;for(e=0;e<d.length;e++){f=d[e][0],g=d[e][1],h=d[e][2];if(g="\
"==null||h===null)continue;var m=g,n=g+h,o=f-(a.gantt.centered?b/2:0),p=f+b-(a.gantt.centered?b/2:0);if(n<"\
"i.min||m>i.max||p<j.min||o>j.max)continue;m<i.min&&(m=i.min),n>i.max&&(n=i.max),o<j.min&&(o=j.min),p>j.ma"\
"x&&(p=j.max);var q=i.d2p(n)-i.d2p(m)-(i.d2p(n)+l<=this.plotWidth?0:l),r=j.d2p(o)-j.d2p(p)-(j.d2p(o)+l<=th"\
"is.plotHeight?0:l);k.fillStyle=\"rgba(0,0,0,0.05)\",k.fillRect(Math.min(i.d2p(m)+l,this.plotWidth),Math.m"\
"in(j.d2p(p)+l,this.plotHeight),q,r)}},extendXRange:function(a){if(a.options.max===null){var b=a.min,c=a.m"\
"ax,d,e,f,g,h,i={},j={},k=null;for(d=0;d<this.series.length;++d){g=this.series[d],h=g.gantt;if(h.show&&g.x"\
"axis==a){for(e=0;e<g.data.length;e++)h.show&&(y=g.data[e][0]+\"\",i[y]=Math.max(i[y]||0,g.data[e][1]+g.da"\
"ta[e][2]),k=g);for(e in i)c=Math.max(i[e],c)}}a.lastSerie=k,a.max=c,a.min=b}},extendYRange:function(a){if"\
"(a.options.max===null){var b=Number.MIN_VALUE,c=Number.MAX_VALUE,d,e,f,g,h={},i={},j=null;for(d=0;d<this."\
"series.length;++d){f=this.series[d],g=f.gantt;if(g.show&&!f.hide&&f.yaxis==a){var k=Number.MIN_VALUE,l=Nu"\
"mber.MAX_VALUE;for(e=0;e<f.data.length;e++)k=Math.max(k,f.data[e][0]),l=Math.min(l,f.data[e][0]);g.center"\
"ed?(b=Math.max(k+.5,b),c=Math.min(l-.5,c)):(b=Math.max(k+1,b),c=Math.min(l,c)),g.barWidth+k>b&&(b=a.max+g"\
".barWidth)}}a.lastSerie=j,a.max=b,a.min=c,a.tickSize=Flotr.getTickSize(a.options.noTicks,c,b,a.options.ti"\
"ckDecimals)}}}),function(){function a(a){return typeof a==\"object\"&&a.constructor&&(Image?!0:a.construc"\
"tor===Image)}Flotr.defaultMarkerFormatter=function(a){return Math.round(a.y*100)/100+\"\"},Flotr.addType("\
"\"markers\",{options:{show:!1,lineWidth:1,color:\"#000000\",fill:!1,fillColor:\"#FFFFFF\",fillOpacity:.4,"\
"stroke:!1,position:\"ct\",verticalMargin:0,labelFormatter:Flotr.defaultMarkerFormatter,fontSize:Flotr.def"\
"aultOptions.fontSize,stacked:!1,stackingType:\"b\",horizontal:!1},stack:{positive:[],negative:[],values:["\
"]},draw:function(a){function m(a,b){return g=d.negative[a]||0,f=d.positive[a]||0,b>0?(d.positive[a]=g+b,g"\
"+b):(d.negative[a]=f+b,f+b)}var b=a.data,c=a.context,d=a.stacked?a.stack:!1,e=a.stackingType,f,g,h,i,j,k,"\
"l;c.save(),c.lineJoin=\"round\",c.lineWidth=a.lineWidth,c.strokeStyle=\"rgba(0,0,0,0.5)\",c.fillStyle=a.f"\
"illStyle;for(i=0;i<b.length;++i)j=b[i][0],k=b[i][1],d&&(e==\"b\"?a.horizontal?k=m(k,j):j=m(j,k):e==\"a\"&"\
"&(h=d.values[j]||0,d.values[j]=h+k,k=h+k)),l=a.labelFormatter({x:j,y:k,index:i,data:b}),this.plot(a.xScal"\
"e(j),a.yScale(k),l,a);c.restore()},plot:function(b,c,d,e){var f=e.context;if(a(d)&&!d.complete)throw\"Mar"\
"ker image not loaded.\";this._plot(b,c,d,e)},_plot:function(b,c,d,e){var f=e.context,g=2,h=b,i=c,j;a(d)?j"\
"={height:d.height,width:d.width}:j=e.text.canvas(d),j.width=Math.floor(j.width+g*2),j.height=Math.floor(j"\
".height+g*2),e.position.indexOf(\"c\")!=-1?h-=j.width/2+g:e.position.indexOf(\"l\")!=-1&&(h-=j.width),e.p"\
"osition.indexOf(\"m\")!=-1?i-=j.height/2+g:e.position.indexOf(\"t\")!=-1?i-=j.height+e.verticalMargin:i+="\
"e.verticalMargin,h=Math.floor(h)+.5,i=Math.floor(i)+.5,e.fill&&f.fillRect(h,i,j.width,j.height),e.stroke&"\
"&f.strokeRect(h,i,j.width,j.height),a(d)?f.drawImage(d,h+g,i+g):Flotr.drawText(f,d,h+g,i+g,{textBaseline:"\
"\"top\",textAlign:\"left\",size:e.fontSize,color:e.color})}})}(),function(){var a=Flotr._;Flotr.defaultPi"\
"eLabelFormatter=function(a,b){return(100*b/a).toFixed(2)+\"%\"},Flotr.addType(\"pie\",{options:{show:!1,l"\
"ineWidth:1,fill:!0,fillColor:null,fillOpacity:.6,explode:6,sizeRatio:.6,startAngle:Math.PI/4,labelFormatt"\
"er:Flotr.defaultPieLabelFormatter,pie3D:!1,pie3DviewAngle:Math.PI/2*.8,pie3DspliceThickness:20},draw:func"\
"tion(a){var b=a.data,c=a.context,d=c.canvas,e=a.lineWidth,f=a.shadowSize,g=a.sizeRatio,h=a.height,i=a.wid"\
"th,j=a.explode,k=a.color,l=a.fill,m=a.fillStyle,n=Math.min(d.width,d.height)*g/2,o=b[0][1],p=[],q=1,r=Mat"\
"h.PI*2*o/this.total,s=this.startAngle||2*Math.PI*a.startAngle,t=s+r,u=s+r/2,v=a.labelFormatter(this.total"\
",o),w=j+n+4,x=Math.cos(u)*w,y=Math.sin(u)*w,z=x<0?\"right\":\"left\",A=y>0?\"top\":\"bottom\",B,C,D,x,y;c"\
".save(),c.translate(i/2,h/2),c.scale(1,q),C=Math.cos(u)*j,D=Math.sin(u)*j,f>0&&(this.plotSlice(C+f,D+f,n,"\
"s,t,c),l&&(c.fillStyle=\"rgba(0,0,0,0.1)\",c.fill())),this.plotSlice(C,D,n,s,t,c),l&&(c.fillStyle=m,c.fil"\
"l()),c.lineWidth=e,c.strokeStyle=k,c.stroke(),B={size:a.fontSize*1.2,color:a.fontColor,weight:1.5},v&&(a."\
"htmlText||!a.textEnabled?(divStyle=\"position:absolute;\"+A+\":\"+(h/2+(A===\"top\"?y:-y))+\"px;\",divSty"\
"le+=z+\":\"+(i/2+(z===\"right\"?-x:x))+\"px;\",p.push('<div style=\"',divStyle,'\" class=\"flotr-grid-lab"\
"el\">',v,\"</div>\")):(B.textAlign=z,B.textBaseline=A,Flotr.drawText(c,v,x,y,B)));if(a.htmlText||!a.textE"\
"nabled){var E=Flotr.DOM.node('<div style=\"color:'+a.fontColor+'\" class=\"flotr-labels\"></div>');Flotr."\
"DOM.insert(E,p.join(\"\")),Flotr.DOM.insert(a.element,E)}c.restore(),this.startAngle=t,this.slices=this.s"\
"lices||[],this.slices.push({radius:Math.min(d.width,d.height)*g/2,x:C,y:D,explode:j,start:s,end:t})},plot"\
"Slice:function(a,b,c,d,e,f){f.beginPath(),f.moveTo(a,b),f.arc(a,b,c,d,e,!1),f.lineTo(a,b),f.closePath()},"\
"hit:function(a){var b=a.data[0],c=a.args,d=a.index,e=c[0],f=c[1],g=this.slices[d],h=e.relX-a.width/2,i=e."\
"relY-a.height/2,j=Math.sqrt(h*h+i*i),k=Math.atan(i/h),l=Math.PI*2,m=g.explode||a.explode,n=g.start%l,o=g."\
"end%l;h<0?k+=Math.PI:h>0&&i<0&&(k+=l),j<g.radius+m&&j>m&&(n>=o&&(k<o||k>n)||k>n&&k<o)&&(f.x=b[0],f.y=b[1]"\
",f.sAngle=n,f.eAngle=o,f.index=0,f.seriesIndex=d,f.fraction=b[1]/this.total)},drawHit:function(a){var b=a"\
".context,c=this.slices[a.args.seriesIndex];b.save(),b.translate(a.width/2,a.height/2),this.plotSlice(c.x,"\
"c.y,c.radius,c.start,c.end,b),b.stroke(),b.restore()},clearHit:function(a){var b=a.context,c=this.slices["\
"a.args.seriesIndex],d=2*a.lineWidth,e=c.radius+d;b.save(),b.translate(a.width/2,a.height/2),b.clearRect(c"\
".x-e,c.y-e,2*e+d,2*e+d),b.restore()},extendYRange:function(a,b){this.total=(this.total||0)+b[0][1]}})}(),"\
"Flotr.addType(\"points\",{options:{show:!1,radius:3,lineWidth:2,fill:!0,fillColor:\"#FFFFFF\",fillOpacity"\
":.4},draw:function(a){var b=a.context,c=a.lineWidth,d=a.shadowSize;b.save(),d>0&&(b.lineWidth=d/2,b.strok"\
"eStyle=\"rgba(0,0,0,0.1)\",this.plot(a,d/2+b.lineWidth/2),b.strokeStyle=\"rgba(0,0,0,0.2)\",this.plot(a,b"\
".lineWidth/2)),b.lineWidth=a.lineWidth,b.strokeStyle=a.color,b.fillStyle=a.fillColor||a.color,this.plot(a"\
"),b.restore()},plot:function(a,b){var c=a.data,d=a.context,e=a.xScale,f=a.yScale,g,h,i;for(g=c.length-1;g"\
">-1;--g){i=c[g][1];if(i===null)continue;h=e(c[g][0]),i=f(i);if(h<0||h>a.width||i<0||i>a.height)continue;d"\
".beginPath(),b?d.arc(h,i+b,a.radius,0,Math.PI,!1):(d.arc(h,i,a.radius,0,2*Math.PI,!0),a.fill&&d.fill()),d"\
".stroke(),d.closePath()}}}),Flotr.addType(\"radar\",{options:{show:!1,lineWidth:2,fill:!0,fillOpacity:.4,"\
"radiusRatio:.9},draw:function(a){var b=a.context,c=a.shadowSize;b.save(),b.translate(a.width/2,a.height/2"\
"),b.lineWidth=a.lineWidth,b.fillStyle=\"rgba(0,0,0,0.05)\",b.strokeStyle=\"rgba(0,0,0,0.05)\",this.plot(a"\
",c/2),b.strokeStyle=\"rgba(0,0,0,0.1)\",this.plot(a,c/4),b.strokeStyle=a.color,b.fillStyle=a.fillStyle,th"\
"is.plot(a),b.restore()},plot:function(a,b){var c=a.data,d=a.context,e=Math.min(a.height,a.width)*a.radius"\
"Ratio/2,f=2*Math.PI/c.length,g=-Math.PI/2,h,i;b=b||0,d.beginPath();for(h=0;h<c.length;++h)i=c[h][1]/this."\
"max,d[h===0?\"moveTo\":\"lineTo\"](Math.cos(h*f+g)*e*i+b,Math.sin(h*f+g)*e*i+b);d.closePath(),a.fill&&d.f"\
"ill(),d.stroke()},extendYRange:function(a,b){this.max=Math.max(a.max,this.max||-Number.MAX_VALUE)}}),Flot"\
"r.addType(\"timeline\",{options:{show:!1,lineWidth:1,barWidth:.2,fill:!0,fillColor:null,fillOpacity:.4,ce"\
"ntered:!0},draw:function(a){var b=a.context;b.save(),b.lineJoin=\"miter\",b.lineWidth=a.lineWidth,b.strok"\
"eStyle=a.color,b.fillStyle=a.fillStyle,this.plot(a),b.restore()},plot:function(a){var b=a.data,c=a.contex"\
"t,d=a.xScale,e=a.yScale,f=a.barWidth,g=a.lineWidth,h;Flotr._.each(b,function(a){var b=a[0],h=a[1],i=a[2],"\
"j=f,k=Math.ceil(d(b)),l=Math.ceil(d(b+i))-k,m=Math.round(e(h)),n=Math.round(e(h-j))-m,o=k-g/2,p=Math.roun"\
"d(m-n/2)-g/2;c.strokeRect(o,p,l,n),c.fillRect(o,p,l,n)})},extendRange:function(a){var b=a.data,c=a.xaxis,"\
"d=a.yaxis,e=a.timeline.barWidth;c.options.min===null&&(c.min=c.datamin-e/2);if(c.options.max===null){var "\
"f=c.max;Flotr._.each(b,function(a){f=Math.max(f,a[0]+a[2])},this),c.max=f+e/2}d.options.min===null&&(d.mi"\
"n=d.datamin-e),d.options.min===null&&(d.max=d.datamax+e)}}),function(){var a=Flotr.DOM;Flotr.addPlugin(\""\
"crosshair\",{options:{mode:null,color:\"#FF0000\",hideCursor:!0},callbacks:{\"flotr:mousemove\":function("\
"a,b){this.options.crosshair.mode&&(this.crosshair.clearCrosshair(),this.crosshair.drawCrosshair(b))}},dra"\
"wCrosshair:function(b){var c=this.octx,d=this.options.crosshair,e=this.plotOffset,f=e.left+b.relX+.5,g=e."\
"top+b.relY+.5;if(b.relX<0||b.relY<0||b.relX>this.plotWidth||b.relY>this.plotHeight){this.el.style.cursor="\
"null,a.removeClass(this.el,\"flotr-crosshair\");return}d.hideCursor&&(this.el.style.cursor=\"none\",a.add"\
"Class(this.el,\"flotr-crosshair\")),c.save(),c.strokeStyle=d.color,c.lineWidth=1,c.beginPath(),d.mode.ind"\
"exOf(\"x\")!=-1&&(c.moveTo(f,e.top),c.lineTo(f,e.top+this.plotHeight)),d.mode.indexOf(\"y\")!=-1&&(c.move"\
"To(e.left,g),c.lineTo(e.left+this.plotWidth,g)),c.stroke(),c.restore()},clearCrosshair:function(){var a=t"\
"his.plotOffset,b=this.lastMousePos,c=this.octx;b&&(c.clearRect(b.relX+a.left,a.top,1,this.plotHeight+1),c"\
".clearRect(a.left,b.relY+a.top,this.plotWidth+1,1))}})}(),function(){function c(a,b,c,d){var e=\"image/\""\
"+a,f=b.toDataURL(e),g=new Image;return g.src=f,g}var a=Flotr.DOM,b=Flotr._;Flotr.addPlugin(\"download\",{"\
"saveImage:function(d,e,f,g){var h=null;if(Flotr.isIE&&Flotr.isIE<9)return h=\"<html><body>\"+this.canvas."\
"firstChild.innerHTML+\"</body></html>\",window.open().document.write(h);if(d!==\"jpeg\"&&d!==\"png\")retu"\
"rn;h=c(d,this.canvas,e,f);if(!b.isElement(h)||!g)return window.open(h.src);this.download.restoreCanvas(),"\
"a.hide(this.canvas),a.hide(this.overlay),a.setStyles({position:\"absolute\"}),a.insert(this.el,h),this.sa"\
"veImageElement=h},restoreCanvas:function(){a.show(this.canvas),a.show(this.overlay),this.saveImageElement"\
"&&this.el.removeChild(this.saveImageElement),this.saveImageElement=null}})}(),function(){var a=Flotr.Even"\
"tAdapter,b=Flotr._;Flotr.addPlugin(\"graphGrid\",{callbacks:{\"flotr:beforedraw\":function(){this.graphGr"\
"id.drawGrid()},\"flotr:afterdraw\":function(){this.graphGrid.drawOutline()}},drawGrid:function(){function"\
" p(a){for(n=0;n<a.length;++n){var b=a[n].v/l.max;for(o=0;o<=u;++o)c[o===0?\"moveTo\":\"lineTo\"](Math.cos"\
"(o*v+w)*t*b,Math.sin(o*v+w)*t*b)}}function q(a,d){b.each(b.pluck(a,\"v\"),function(a){if(a<=l.min||a>=l.m"\
"ax||(a==l.min||a==l.max)&&e.outlineWidth)return;d(Math.floor(l.d2p(a))+c.lineWidth/2)})}function r(a){c.m"\
"oveTo(a,0),c.lineTo(a,j)}function s(a){c.moveTo(0,a),c.lineTo(k,a)}var c=this.ctx,d=this.options,e=d.grid"\
",f=e.verticalLines,g=e.horizontalLines,h=e.minorVerticalLines,i=e.minorHorizontalLines,j=this.plotHeight,"\
"k=this.plotWidth,l,m,n,o;(f||h||g||i)&&a.fire(this.el,\"flotr:beforegrid\",[this.axes.x,this.axes.y,d,thi"\
"s]),c.save(),c.lineWidth=1,c.strokeStyle=e.tickColor;if(e.circular){c.translate(this.plotOffset.left+k/2,"\
"this.plotOffset.top+j/2);var t=Math.min(j,k)*d.radar.radiusRatio/2,u=this.axes.x.ticks.length,v=2*(Math.P"\
"I/u),w=-Math.PI/2;c.beginPath(),l=this.axes.y,g&&p(l.ticks),i&&p(l.minorTicks),f&&b.times(u,function(a){c"\
".moveTo(0,0),c.lineTo(Math.cos(a*v+w)*t,Math.sin(a*v+w)*t)}),c.stroke()}else c.translate(this.plotOffset."\
"left,this.plotOffset.top),e.backgroundColor&&(c.fillStyle=this.processColor(e.backgroundColor,{x1:0,y1:0,"\
"x2:k,y2:j}),c.fillRect(0,0,k,j)),c.beginPath(),l=this.axes.x,f&&q(l.ticks,r),h&&q(l.minorTicks,r),l=this."\
"axes.y,g&&q(l.ticks,s),i&&q(l.minorTicks,s),c.stroke();c.restore(),(f||h||g||i)&&a.fire(this.el,\"flotr:a"\
"ftergrid\",[this.axes.x,this.axes.y,d,this])},drawOutline:function(){var a=this,b=a.options,c=b.grid,d=c."\
"outline,e=a.ctx,f=c.backgroundImage,g=a.plotOffset,h=g.left,j=g.top,k=a.plotWidth,l=a.plotHeight,m,n,o,p,"\
"q,r;if(!c.outlineWidth)return;e.save();if(c.circular){e.translate(h+k/2,j+l/2);var s=Math.min(l,k)*b.rada"\
"r.radiusRatio/2,t=this.axes.x.ticks.length,u=2*(Math.PI/t),v=-Math.PI/2;e.beginPath(),e.lineWidth=c.outli"\
"neWidth,e.strokeStyle=c.color,e.lineJoin=\"round\";for(i=0;i<=t;++i)e[i===0?\"moveTo\":\"lineTo\"](Math.c"\
"os(i*u+v)*s,Math.sin(i*u+v)*s);e.stroke()}else{e.translate(h,j);var w=c.outlineWidth,x=.5-w+(w+1)%2/2,y=\"\
""lineTo\",z=\"moveTo\";e.lineWidth=w,e.strokeStyle=c.color,e.lineJoin=\"miter\",e.beginPath(),e.moveTo(x,"\
"x),k-=w/2%1,l+=w/2,e[d.indexOf(\"n\")!==-1?y:z](k,x),e[d.indexOf(\"e\")!==-1?y:z](k,l),e[d.indexOf(\"s\")"\
"!==-1?y:z](x,l),e[d.indexOf(\"w\")!==-1?y:z](x,x),e.stroke(),e.closePath()}e.restore(),f&&(o=f.src||f,p=("\
"parseInt(f.left,10)||0)+g.left,q=(parseInt(f.top,10)||0)+g.top,n=new Image,n.onload=function(){e.save(),f"\
".alpha&&(e.globalAlpha=f.alpha),e.globalCompositeOperation=\"destination-over\",e.drawImage(n,0,0,n.width"\
",n.height,p,q,k,l),e.restore()},n.src=o)}})}(),function(){var a=Flotr.DOM,b=Flotr._,c=Flotr,d=\"opacity:0"\
".7;background-color:#000;color:#fff;display:none;position:absolute;padding:2px 8px;-moz-border-radius:4px"\
";border-radius:4px;white-space:nowrap;\";Flotr.addPlugin(\"hit\",{callbacks:{\"flotr:mousemove\":function"\
"(a,b){this.hit.track(b)},\"flotr:click\":function(a){this.hit.track(a)},\"flotr:mouseout\":function(){thi"\
"s.hit.clearHit()}},track:function(a){(this.options.mouse.track||b.any(this.series,function(a){return a.mo"\
"use&&a.mouse.track}))&&this.hit.hit(a)},executeOnType:function(a,d,e){function h(a,h){b.each(b.keys(c.gra"\
"phTypes),function(b){a[b]&&a[b].show&&this[b][d]&&(g=this.getOptions(a,b),g.fill=!!a.mouse.fillColor,g.fi"\
"llStyle=this.processColor(a.mouse.fillColor||\"#ffffff\",{opacity:a.mouse.fillOpacity}),g.color=a.mouse.l"\
"ineColor,g.context=this.octx,g.index=h,e&&(g.args=e),this[b][d].call(this[b],g),f=!0)},this)}var f=!1,g;r"\
"eturn b.isArray(a)||(a=[a]),b.each(a,h,this),f},drawHit:function(a){var b=this.octx,c=a.series;if(c.mouse"\
".lineColor){b.save(),b.lineWidth=c.points?c.points.lineWidth:1,b.strokeStyle=c.mouse.lineColor,b.fillStyl"\
"e=this.processColor(c.mouse.fillColor||\"#ffffff\",{opacity:c.mouse.fillOpacity}),b.translate(this.plotOf"\
"fset.left,this.plotOffset.top);if(!this.hit.executeOnType(c,\"drawHit\",a)){var d=a.xaxis,e=a.yaxis;b.beg"\
"inPath(),b.arc(d.d2p(a.x),e.d2p(a.y),c.points.radius||c.mouse.radius,0,2*Math.PI,!0),b.fill(),b.stroke(),"\
"b.closePath()}b.restore(),this.clip(b)}this.prevHit=a},clearHit:function(){var b=this.prevHit,c=this.octx"\
",d=this.plotOffset;c.save(),c.translate(d.left,d.top);if(b){if(!this.hit.executeOnType(b.series,\"clearHi"\
"t\",this.prevHit)){var e=b.series,f=e.points?e.points.lineWidth:1;offset=(e.points.radius||e.mouse.radius"\
")+f,c.clearRect(b.xaxis.d2p(b.x)-offset,b.yaxis.d2p(b.y)-offset,offset*2,offset*2)}a.hide(this.mouseTrack"\
"),this.prevHit=null}c.restore()},hit:function(a){var c=this.options,d=this.prevHit,e,f,g,h,i,j,k,l;if(thi"\
"s.series.length===0)return;n={relX:a.relX,relY:a.relY,absX:a.absX,absY:a.absY};if(c.mouse.trackY&&!c.mous"\
"e.trackAll&&this.hit.executeOnType(this.series,\"hit\",[a,n]))b.isUndefined(n.seriesIndex)||(i=this.serie"\
"s[n.seriesIndex],n.series=i,n.mouse=i.mouse,n.xaxis=i.xaxis,n.yaxis=i.yaxis);else{e=this.hit.closest(a);i"\
"f(e){e=c.mouse.trackY?e.point:e.x,h=e.seriesIndex,i=this.series[h],k=i.xaxis,l=i.yaxis,f=2*i.mouse.sensib"\
"ility;if(c.mouse.trackAll||e.distanceX<f/k.scale&&(!c.mouse.trackY||e.distanceY<f/l.scale))n.series=i,n.x"\
"axis=i.xaxis,n.yaxis=i.yaxis,n.mouse=i.mouse,n.x=e.x,n.y=e.y,n.dist=e.distance,n.index=e.dataIndex,n.seri"\
"esIndex=h}}if(!d||d.index!==n.index||d.seriesIndex!==n.seriesIndex)this.hit.clearHit(),n.series&&n.mouse&"\
"&n.mouse.track&&(this.hit.drawMouseTrack(n),this.hit.drawHit(n),Flotr.EventAdapter.fire(this.el,\"flotr:h"\
"it\",[n,this]))},closest:function(a){function v(a){a.distance=m,a.distanceX=n,a.distanceY=o,a.seriesIndex"\
"=t,a.dataIndex=u,a.x=r,a.y=s}var b=this.series,c=this.options,d=a.relX,e=a.relY,f=Number.MAX_VALUE,g=Numb"\
"er.MAX_VALUE,h={},i={},j=!1,k,l,m,n,o,p,q,r,s,t,u;for(t=0;t<b.length;t++){k=b[t],l=k.data,p=k.xaxis.p2d(d"\
"),q=k.yaxis.p2d(e),l.length&&(j=!0);for(u=l.length;u--;){r=l[u][0],s=l[u][1];if(r===null||s===null)contin"\
"ue;if(r<k.xaxis.min||r>k.xaxis.max)continue;n=Math.abs(r-p),o=Math.abs(s-q),m=n*n+o*o,m<f&&(f=m,v(h)),n<g"\
"&&(g=n,v(i))}}return j?{point:h,x:i}:!1},drawMouseTrack:function(b){var c=\"\",e=b.series,f=b.mouse.posit"\
"ion,g=b.mouse.margin,h=d,i=this.mouseTrack,j=this.plotOffset,k=j.left,l=j.right,m=j.bottom,n=j.top,o=b.mo"\
"use.trackDecimals,p=this.options;i||(i=a.node('<div class=\"flotr-mouse-value\"></div>'),this.mouseTrack="\
"i,a.insert(this.el,i));if(!b.mouse.relative)f.charAt(0)==\"n\"?c+=\"top:\"+(g+n)+\"px;bottom:auto;\":f.ch"\
"arAt(0)==\"s\"&&(c+=\"bottom:\"+(g+m)+\"px;top:auto;\"),f.charAt(1)==\"e\"?c+=\"right:\"+(g+l)+\"px;left:"\
"auto;\":f.charAt(1)==\"w\"&&(c+=\"left:\"+(g+k)+\"px;right:auto;\");else if(e.bars.show)c+=\"bottom:\"+(g"\
"-n-b.yaxis.d2p(b.y/2)+this.canvasHeight)+\"px;top:auto;\",c+=\"left:\"+(g+k+b.xaxis.d2p(b.x-p.bars.barWid"\
"th/2))+\"px;right:auto;\";else if(e.pie.show){var q={x:this.plotWidth/2,y:this.plotHeight/2},r=Math.min(t"\
"his.canvasWidth,this.canvasHeight)*e.pie.sizeRatio/2,s=b.sAngle<b.eAngle?(b.sAngle+b.eAngle)/2:(b.sAngle+"\
"b.eAngle+2*Math.PI)/2;c+=\"bottom:\"+(g-n-q.y-Math.sin(s)*r/2+this.canvasHeight)+\"px;top:auto;\",c+=\"le"\
"ft:\"+(g+k+q.x+Math.cos(s)*r/2)+\"px;right:auto;\"}else f.charAt(0)==\"n\"?c+=\"bottom:\"+(g-n-b.yaxis.d2"\
"p(b.y)+this.canvasHeight)+\"px;top:auto;\":f.charAt(0)==\"s\"&&(c+=\"top:\"+(g+n+b.yaxis.d2p(b.y))+\"px;b"\
"ottom:auto;\"),f.charAt(1)==\"e\"?c+=\"left:\"+(g+k+b.xaxis.d2p(b.x))+\"px;right:auto;\":f.charAt(1)==\"w"\
"\"&&(c+=\"right:\"+(g-k-b.xaxis.d2p(b.x)+this.canvasWidth)+\"px;left:auto;\");h+=c,i.style.cssText=h;if(!"\
"o||o<0)o=0;i.innerHTML=b.mouse.trackFormatter({x:b.x.toFixed(o),y:b.y.toFixed(o),series:b.series,index:b."\
"index,nearest:b,fraction:b.fraction}),a.show(i)}})}(),function(){function a(a,b){return a.which?a.which=="\
"=1:a.button===0||a.button===1}function b(a,b){return Math.min(Math.max(0,a),b.plotWidth-1)}function c(a,b"\
"){return Math.min(Math.max(0,a),b.plotHeight)}var d=Flotr.DOM,e=Flotr.EventAdapter,f=Flotr._;Flotr.addPlu"\
"gin(\"selection\",{options:{pinchOnly:null,mode:null,color:\"#B6D9FF\",fps:20},callbacks:{\"flotr:mouseup"\
"\":function(a){var b=this.options.selection,c=this.selection,d=this.getEventPosition(a);if(!b||!b.mode)re"\
"turn;c.interval&&clearInterval(c.interval),this.multitouches?c.updateSelection():b.pinchOnly||c.setSelect"\
"ionPos(c.selection.second,d),c.clearSelection(),c.selecting&&c.selectionIsSane()&&(c.drawSelection(),c.fi"\
"reSelectEvent(),this.ignoreClick=!0)},\"flotr:mousedown\":function(b){var c=this.options.selection,d=this"\
".selection,e=this.getEventPosition(b);if(!c||!c.mode)return;if(!c.mode||!a(b)&&f.isUndefined(b.touches))r"\
"eturn;c.pinchOnly||d.setSelectionPos(d.selection.first,e),d.interval&&clearInterval(d.interval),this.last"\
"MousePos.pageX=null,d.selecting=!1,d.interval=setInterval(f.bind(d.updateSelection,this),1e3/c.fps)},\"fl"\
"otr:destroy\":function(a){clearInterval(this.selection.interval)}},getArea:function(){var a=this.selectio"\
"n.selection,b=a.first,c=a.second;return{x1:Math.min(b.x,c.x),x2:Math.max(b.x,c.x),y1:Math.min(b.y,c.y),y2"\
":Math.max(b.y,c.y)}},selection:{first:{x:-1,y:-1},second:{x:-1,y:-1}},prevSelection:null,interval:null,fi"\
"reSelectEvent:function(a){var b=this.axes,c=this.selection.selection,d,f,g,h;a=a||\"select\",d=b.x.p2d(c."\
"first.x),f=b.x.p2d(c.second.x),g=b.y.p2d(c.first.y),h=b.y.p2d(c.second.y),e.fire(this.el,\"flotr:\"+a,[{x"\
"1:Math.min(d,f),y1:Math.min(g,h),x2:Math.max(d,f),y2:Math.max(g,h),xfirst:d,xsecond:f,yfirst:g,ysecond:h}"\
",this])},setSelection:function(a,d){var e=this.options,f=this.axes.x,g=this.axes.y,h=g.scale,i=f.scale,j="\
"e.selection.mode.indexOf(\"x\")!=-1,k=e.selection.mode.indexOf(\"y\")!=-1,l=this.selection.selection;this"\
".selection.clearSelection(),l.first.y=c(j&&!k?0:(g.max-a.y1)*h,this),l.second.y=c(j&&!k?this.plotHeight-1"\
":(g.max-a.y2)*h,this),l.first.x=b(k&&!j?0:a.x1,this),l.second.x=b(k&&!j?this.plotWidth:a.x2,this),this.se"\
"lection.drawSelection(),d||this.selection.fireSelectEvent()},setSelectionPos:function(a,d){var e=this.opt"\
"ions.selection.mode,f=this.selection.selection;e.indexOf(\"x\")==-1?a.x=a==f.first?0:this.plotWidth:a.x=b"\
"(d.relX,this),e.indexOf(\"y\")==-1?a.y=a==f.first?0:this.plotHeight-1:a.y=c(d.relY,this)},drawSelection:f"\
"unction(){this.selection.fireSelectEvent(\"selecting\");var a=this.selection.selection,b=this.octx,c=this"\
".options,d=this.plotOffset,e=this.selection.prevSelection;if(e&&a.first.x==e.first.x&&a.first.y==e.first."\
"y&&a.second.x==e.second.x&&a.second.y==e.second.y)return;b.save(),b.strokeStyle=this.processColor(c.selec"\
"tion.color,{opacity:.8}),b.lineWidth=1,b.lineJoin=\"miter\",b.fillStyle=this.processColor(c.selection.col"\
"or,{opacity:.4}),this.selection.prevSelection={first:{x:a.first.x,y:a.first.y},second:{x:a.second.x,y:a.s"\
"econd.y}};var f=Math.min(a.first.x,a.second.x),g=Math.min(a.first.y,a.second.y),h=Math.abs(a.second.x-a.f"\
"irst.x),i=Math.abs(a.second.y-a.first.y);b.fillRect(f+d.left+.5,g+d.top+.5,h,i),b.strokeRect(f+d.left+.5,"\
"g+d.top+.5,h,i),b.restore()},updateSelection:function(){if(!this.lastMousePos.pageX)return;this.selection"\
".selecting=!0;if(this.multitouches)this.selection.setSelectionPos(this.selection.selection.first,this.get"\
"EventPosition(this.multitouches[0])),this.selection.setSelectionPos(this.selection.selection.second,this."\
"getEventPosition(this.multitouches[1]));else{if(this.options.selection.pinchOnly)return;this.selection.se"\
"tSelectionPos(this.selection.selection.second,this.lastMousePos)}this.selection.clearSelection(),this.sel"\
"ection.selectionIsSane()&&this.selection.drawSelection()},clearSelection:function(){if(!this.selection.pr"\
"evSelection)return;var a=this.selection.prevSelection,b=1,c=this.plotOffset,d=Math.min(a.first.x,a.second"\
".x),e=Math.min(a.first.y,a.second.y),f=Math.abs(a.second.x-a.first.x),g=Math.abs(a.second.y-a.first.y);th"\
"is.octx.clearRect(d+c.left-b+.5,e+c.top-b,f+2*b+.5,g+2*b+.5),this.selection.prevSelection=null},selection"\
"IsSane:function(){var a=this.selection.selection;return Math.abs(a.second.x-a.first.x)>=5||Math.abs(a.sec"\
"ond.y-a.first.y)>=5}})}(),function(){var a=Flotr.DOM;Flotr.addPlugin(\"labels\",{callbacks:{\"flotr:after"\
"draw\":function(){this.labels.draw()}},draw:function(){function s(a,b,d){var e=d?b.minorTicks:b.ticks,f=b"\
".orientation===1,h=b.n===1,k,m;k={color:b.options.color||o.grid.color,angle:Flotr.toRad(b.options.labelsA"\
"ngle),textBaseline:\"middle\"};for(l=0;l<e.length&&(d?b.options.showMinorLabels:b.options.showLabels);++l"\
"){c=e[l],c.label+=\"\";if(!c.label||!c.label.length)continue;x=Math.cos(l*i+j)*g,y=Math.sin(l*i+j)*g,k.te"\
"xtAlign=f?Math.abs(x)<.1?\"center\":x<0?\"right\":\"left\":\"left\",Flotr.drawText(p,c.label,f?x:3,f?y:-("\
"b.ticks[l].v/b.max)*(g-o.fontSize),k)}}function t(a,b,d,e){function j(a){return a.options.showLabels&&a.u"\
"sed}function k(a,b,c,d){return a.plotOffset.left+(b?d:c?-o.grid.labelMargin:o.grid.labelMargin+a.plotWidt"\
"h)}function m(a,b,c,d){return a.plotOffset.top+(b?o.grid.labelMargin:d)+(b&&c?a.plotHeight:0)}var f=b.ori"\
"entation===1,g=b.n===1,h,i;h={color:b.options.color||o.grid.color,textAlign:d,textBaseline:e,angle:Flotr."\
"toRad(b.options.labelsAngle)},h=Flotr.getBestTextAlign(h.angle,h);for(l=0;l<b.ticks.length&&j(b);++l){c=b"\
".ticks[l];if(!c.label||!c.label.length)continue;i=b.d2p(c.v);if(i<0||i>(f?a.plotWidth:a.plotHeight))conti"\
"nue;Flotr.drawText(p,c.label,k(a,f,g,i),m(a,f,g,i),h),!f&&!g&&(p.save(),p.strokeStyle=h.color,p.beginPath"\
"(),p.moveTo(a.plotOffset.left+a.plotWidth-8,a.plotOffset.top+b.d2p(c.v)),p.lineTo(a.plotOffset.left+a.plo"\
"tWidth,a.plotOffset.top+b.d2p(c.v)),p.stroke(),p.restore())}}function u(a,b){var d=b.orientation===1,e=b."\
"n===1,g=\"\",h,i,j,k=a.plotOffset;!d&&!e&&(p.save(),p.strokeStyle=b.options.color||o.grid.color,p.beginPa"\
"th());if(b.options.showLabels&&(e?!0:b.used))for(l=0;l<b.ticks.length;++l){c=b.ticks[l];if(!c.label||!c.l"\
"abel.length||(d?k.left:k.top)+b.d2p(c.v)<0||(d?k.left:k.top)+b.d2p(c.v)>(d?a.canvasWidth:a.canvasHeight))"\
"continue;j=k.top+(d?(e?1:-1)*(a.plotHeight+o.grid.labelMargin):b.d2p(c.v)-b.maxLabel.height/2),h=d?k.left"\
"+b.d2p(c.v)-f/2:0,g=\"\",l===0?g=\" first\":l===b.ticks.length-1&&(g=\" last\"),g+=d?\" flotr-grid-label-"\
"x\":\" flotr-grid-label-y\",m+=['<div style=\"position:absolute; text-align:'+(d?\"center\":\"right\")+\""\
"; \",\"top:\"+j+\"px; \",(!d&&!e?\"right:\":\"left:\")+h+\"px; \",\"width:\"+(d?f:(e?k.left:k.right)-o.gr"\
"id.labelMargin)+\"px; \",b.options.color?\"color:\"+b.options.color+\"; \":\" \",'\" class=\"flotr-grid-l"\
"abel'+g+'\">'+c.label+\"</div>\"].join(\" \"),!d&&!e&&(p.moveTo(k.left+a.plotWidth-8,k.top+b.d2p(c.v)),p."\
"lineTo(k.left+a.plotWidth,k.top+b.d2p(c.v)))}}var b,c,d,e,f,g,h,i,j,k,l,m=\"\",n=0,o=this.options,p=this."\
"ctx,q=this.axes,r={size:o.fontSize};for(l=0;l<q.x.ticks.length;++l)q.x.ticks[l].label&&++n;f=this.plotWid"\
"th/n,o.grid.circular&&(p.save(),p.translate(this.plotOffset.left+this.plotWidth/2,this.plotOffset.top+thi"\
"s.plotHeight/2),g=this.plotHeight*o.radar.radiusRatio/2+o.fontSize,h=this.axes.x.ticks.length,i=2*(Math.P"\
"I/h),j=-Math.PI/2,s(this,q.x,!1),s(this,q.x,!0),s(this,q.y,!1),s(this,q.y,!0),p.restore()),!o.HtmlText&&t"\
"his.textEnabled?(t(this,q.x,\"center\",\"top\"),t(this,q.x2,\"center\",\"bottom\"),t(this,q.y,\"right\",\"\
""middle\"),t(this,q.y2,\"left\",\"middle\")):(q.x.options.showLabels||q.x2.options.showLabels||q.y.option"\
"s.showLabels||q.y2.options.showLabels)&&!o.grid.circular&&(m=\"\",u(this,q.x),u(this,q.x2),u(this,q.y),u("\
"this,q.y2),p.stroke(),p.restore(),k=a.create(\"div\"),a.setStyles(k,{fontSize:\"smaller\",color:o.grid.co"\
"lor}),k.className=\"flotr-labels\",a.insert(this.el,k),a.insert(k,m))}})}(),function(){var a=Flotr.DOM,b="\
"Flotr._;Flotr.addPlugin(\"legend\",{options:{show:!0,noColumns:1,labelFormatter:function(a){return a},lab"\
"elBoxBorderColor:\"#CCCCCC\",labelBoxWidth:14,labelBoxHeight:10,labelBoxMargin:5,labelBoxOpacity:.4,conta"\
"iner:null,position:\"nw\",margin:5,backgroundColor:null,backgroundOpacity:.85},callbacks:{\"flotr:afterin"\
"it\":function(){this.legend.insertLegend()}},insertLegend:function(){if(!this.options.legend.show)return;"\
"var c=this.series,d=this.plotOffset,e=this.options,f=e.legend,g=[],h=!1,i=this.ctx,j=b.filter(c,function("\
"a){return a.label&&!a.hide}).length,k=f.position,l=f.margin,m,n,o;if(j)if(!e.HtmlText&&this.textEnabled&&"\
"!f.container){var p={size:e.fontSize*1.1,color:e.grid.color},q=f.labelBoxWidth,r=f.labelBoxHeight,s=f.lab"\
"elBoxMargin,t=d.left+l,u=d.top+l,v=0;for(m=c.length-1;m>-1;--m){if(!c[m].label||c[m].hide)continue;n=f.la"\
"belFormatter(c[m].label),v=Math.max(v,this._text.measureText(n,p).width)}var w=Math.round(q+s*3+v),x=Math"\
".round(j*(s+r)+s);k.charAt(0)==\"s\"&&(u=d.top+this.plotHeight-(l+x)),k.charAt(1)==\"e\"&&(t=d.left+this."\
"plotWidth-(l+w)),o=this.processColor(f.backgroundColor||\"rgb(240,240,240)\",{opacity:f.backgroundOpacity"\
"||.1}),i.fillStyle=o,i.fillRect(t,u,w,x),i.strokeStyle=f.labelBoxBorderColor,i.strokeRect(Flotr.toPixel(t"\
"),Flotr.toPixel(u),w,x);var y=t+s,z=u+s;for(m=0;m<c.length;m++){if(!c[m].label||c[m].hide)continue;n=f.la"\
"belFormatter(c[m].label),i.fillStyle=c[m].color,i.fillRect(y,z,q-1,r-1),i.strokeStyle=f.labelBoxBorderCol"\
"or,i.lineWidth=1,i.strokeRect(Math.ceil(y)-1.5,Math.ceil(z)-1.5,q+2,r+2),Flotr.drawText(i,n,y+q+s,z+r,p),"\
"z+=r+s}}else{for(m=0;m<c.length;++m){if(!c[m].label||c[m].hide)continue;m%f.noColumns===0&&(g.push(h?\"</"\
"tr><tr>\":\"<tr>\"),h=!0);var A=c[m],B=f.labelBoxWidth,C=f.labelBoxHeight,E=A.bars?A.bars.fillOpacity:f.l"\
"abelBoxOpacity,F=\"opacity:\"+E+\";filter:alpha(opacity=\"+E*100+\");\";n=f.labelFormatter(A.label),o=\"b"\
"ackground-color:\"+(A.bars&&A.bars.show&&A.bars.fillColor&&A.bars.fill?A.bars.fillColor:A.color)+\";\",g."\
"push('<td class=\"flotr-legend-color-box\">','<div style=\"border:1px solid ',f.labelBoxBorderColor,';pad"\
"ding:1px\">','<div style=\"width:',B-1,\"px;height:\",C-1,\"px;border:1px solid \",c[m].color,'\">','<div"\
" style=\"width:',B,\"px;height:\",C,\"px;\",\"opacity:.4;\",o,'\"></div>',\"</div>\",\"</div>\",\"</td>\""\
",'<td class=\"flotr-legend-label\">',n,\"</td>\")}h&&g.push(\"</tr>\");if(g.length>0){var G='<table style"\
"=\"font-size:smaller;color:'+e.grid.color+'\">'+g.join(\"\")+\"</table>\";if(f.container)a.insert(f.conta"\
"iner,G);else{var H={position:\"absolute\",\"z-index\":2};k.charAt(0)==\"n\"?(H.top=l+d.top+\"px\",H.botto"\
"m=\"auto\"):k.charAt(0)==\"s\"&&(H.bottom=l+d.bottom+\"px\",H.top=\"auto\"),k.charAt(1)==\"e\"?(H.right=l"\
"+d.right+\"px\",H.left=\"auto\"):k.charAt(1)==\"w\"&&(H.left=l+d.left+\"px\",H.right=\"auto\");var I=a.cr"\
"eate(\"div\"),J;I.className=\"flotr-legend\",a.setStyles(I,H),a.insert(I,G),a.insert(this.el,I);if(!f.bac"\
"kgroundOpacity)return;var K=f.backgroundColor||e.grid.backgroundColor||\"#ffffff\";b.extend(H,a.size(I),{"\
"backgroundColor:K,\"z-index\":1}),H.width+=\"px\",H.height+=\"px\",I=a.create(\"div\"),I.className=\"flot"\
"r-legend-bg\",a.setStyles(I,H),a.opacity(I,f.backgroundOpacity),a.insert(I,\" \"),a.insert(this.el,I)}}}}"\
"})}(),function(){function a(a){if(this.options.spreadsheet.tickFormatter)return this.options.spreadsheet."\
"tickFormatter(a);var b=c.find(this.axes.x.ticks,function(b){return b.v==a});return b?b.label:a}var b=Flot"\
"r.DOM,c=Flotr._;Flotr.addPlugin(\"spreadsheet\",{options:{show:!1,tabGraphLabel:\"Graph\",tabDataLabel:\""\
"Data\",toolbarDownload:\"Download CSV\",toolbarSelectAll:\"Select all\",csvFileSeparator:\",\",decimalSep"\
"arator:\".\",tickFormatter:null,initialTab:\"graph\"},callbacks:{\"flotr:afterconstruct\":function(){if(!"\
"this.options.spreadsheet.show)return;var a=this.spreadsheet,c=b.node('<div class=\"flotr-tabs-group\" sty"\
"le=\"position:absolute;left:0px;width:'+this.canvasWidth+'px\"></div>'),d=b.node('<div style=\"float:left"\
"\" class=\"flotr-tab selected\">'+this.options.spreadsheet.tabGraphLabel+\"</div>\"),e=b.node('<div style"\
"=\"float:left\" class=\"flotr-tab\">'+this.options.spreadsheet.tabDataLabel+\"</div>\"),f;a.tabsContainer"\
"=c,a.tabs={graph:d,data:e},b.insert(c,d),b.insert(c,e),b.insert(this.el,c),f=b.size(e).height+2,this.plot"\
"Offset.bottom+=f,b.setStyles(c,{top:this.canvasHeight-f+\"px\"}),this.observe(d,\"click\",function(){a.sh"\
"owTab(\"graph\")}).observe(e,\"click\",function(){a.showTab(\"data\")}),this.options.spreadsheet.initialT"\
"ab!==\"graph\"&&a.showTab(this.options.spreadsheet.initialTab)}},loadDataGrid:function(){if(this.seriesDa"\
"ta)return this.seriesData;var a=this.series,b={};return c.each(a,function(a,d){c.each(a.data,function(a){"\
"var c=a[0],e=a[1],f=b[c];if(f)f[d+1]=e;else{var g=[];g[0]=c,g[d+1]=e,b[c]=g}})}),this.seriesData=c.sortBy"\
"(b,function(a,b){return parseInt(b,10)}),this.seriesData},constructDataGrid:function(){if(this.spreadshee"\
"t.datagrid)return this.spreadsheet.datagrid;var d=this.series,e=this.spreadsheet.loadDataGrid(),f=[\"<col"\
"group><col />\"],g,h,i,j=['<table class=\"flotr-datagrid\"><tr class=\"first-row\">'];j.push(\"<th>&nbsp;"\
"</th>\"),c.each(d,function(a,b){j.push('<th scope=\"col\">'+(a.label||String.fromCharCode(65+b))+\"</th>\"\
""),f.push(\"<col />\")}),j.push(\"</tr>\"),c.each(e,function(b){j.push(\"<tr>\"),c.times(d.length+1,funct"\
"ion(d){var e=\"td\",f=b[d],g=c.isUndefined(f)?\"\":Math.round(f*1e5)/1e5;if(d===0){e=\"th\";var h=a.call("\
"this,g);h&&(g=h)}j.push(\"<\"+e+(e==\"th\"?' scope=\"row\"':\"\")+\">\"+g+\"</\"+e+\">\")},this),j.push(\"\
""</tr>\")},this),f.push(\"</colgroup>\"),i=b.node(j.join(\"\")),g=b.node('<button type=\"button\" class=\"\
""flotr-datagrid-toolbar-button\">'+this.options.spreadsheet.toolbarDownload+\"</button>\"),h=b.node('<but"\
"ton type=\"button\" class=\"flotr-datagrid-toolbar-button\">'+this.options.spreadsheet.toolbarSelectAll+\"\
""</button>\"),this.observe(g,\"click\",c.bind(this.spreadsheet.downloadCSV,this)).observe(h,\"click\",c.b"\
"ind(this.spreadsheet.selectAllData,this));var k=b.node('<div class=\"flotr-datagrid-toolbar\"></div>');b."\
"insert(k,g),b.insert(k,h);var l=this.canvasHeight-b.size(this.spreadsheet.tabsContainer).height-2,m=b.nod"\
"e('<div class=\"flotr-datagrid-container\" style=\"position:absolute;left:0px;top:0px;width:'+this.canvas"\
"Width+\"px;height:\"+l+'px;overflow:auto;z-index:10\"></div>');return b.insert(m,k),b.insert(m,i),b.inser"\
"t(this.el,m),this.spreadsheet.datagrid=i,this.spreadsheet.container=m,i},showTab:function(a){if(this.spre"\
"adsheet.activeTab===a)return;switch(a){case\"graph\":b.hide(this.spreadsheet.container),b.removeClass(thi"\
"s.spreadsheet.tabs.data,\"selected\"),b.addClass(this.spreadsheet.tabs.graph,\"selected\");break;case\"da"\
"ta\":this.spreadsheet.datagrid||this.spreadsheet.constructDataGrid(),b.show(this.spreadsheet.container),b"\
".addClass(this.spreadsheet.tabs.data,\"selected\"),b.removeClass(this.spreadsheet.tabs.graph,\"selected\""\
");break;default:throw\"Illegal tab name: \"+a}this.spreadsheet.activeTab=a},selectAllData:function(){if(t"\
"his.spreadsheet.tabs){var a,b,c,d,e=this.spreadsheet.constructDataGrid();return this.spreadsheet.showTab("\
"\"data\"),setTimeout(function(){(c=e.ownerDocument)&&(d=c.defaultView)&&d.getSelection&&c.createRange&&(a"\
"=window.getSelection())&&a.removeAllRanges?(b=c.createRange(),b.selectNode(e),a.removeAllRanges(),a.addRa"\
"nge(b)):document.body&&document.body.createTextRange&&(b=document.body.createTextRange())&&(b.moveToEleme"\
"ntText(e),b.select())},0),!0}return!1},downloadCSV:function(){var b=\"\",d=this.series,e=this.options,f=t"\
"his.spreadsheet.loadDataGrid(),g=encodeURIComponent(e.spreadsheet.csvFileSeparator);if(e.spreadsheet.deci"\
"malSeparator===e.spreadsheet.csvFileSeparator)throw\"The decimal separator is the same as the column sepa"\
"rator (\"+e.spreadsheet.decimalSeparator+\")\";c.each(d,function(a,c){b+=g+'\"'+(a.label||String.fromChar"\
"Code(65+c)).replace(/\\\"/g,'\\\\\"')+'\"'}),b+=\"%0D%0A\",b+=c.reduce(f,function(b,c){var d=a.call(this,"\
"c[0])||\"\";d='\"'+(d+\"\").replace(/\\\"/g,'\\\\\"')+'\"';var f=c.slice(1).join(g);return e.spreadsheet."\
"decimalSeparator!==\".\"&&(f=f.replace(/\\./g,e.spreadsheet.decimalSeparator)),b+d+g+f+\"%0D%0A\"},\"\",t"\
"his),Flotr.isIE&&Flotr.isIE<9?(b=b.replace(new RegExp(g,\"g\"),decodeURIComponent(g)).replace(/%0A/g,\"\\"\
"n\").replace(/%0D/g,\"\\r\"),window.open().document.write(b)):window.open(\"data:text/csv,\"+b)}})}(),fun"\
"ction(){var a=Flotr.DOM;Flotr.addPlugin(\"titles\",{callbacks:{\"flotr:afterdraw\":function(){this.titles"\
".drawTitles()}},drawTitles:function(){var b,c=this.options,d=c.grid.labelMargin,e=this.ctx,f=this.axes;if"\
"(!c.HtmlText&&this.textEnabled){var g={size:c.fontSize,color:c.grid.color,textAlign:\"center\"};c.subtitl"\
"e&&Flotr.drawText(e,c.subtitle,this.plotOffset.left+this.plotWidth/2,this.titleHeight+this.subtitleHeight"\
"-2,g),g.weight=1.5,g.size*=1.5,c.title&&Flotr.drawText(e,c.title,this.plotOffset.left+this.plotWidth/2,th"\
"is.titleHeight-2,g),g.weight=1.8,g.size*=.8,f.x.options.title&&f.x.used&&(g.textAlign=f.x.options.titleAl"\
"ign||\"center\",g.textBaseline=\"top\",g.angle=Flotr.toRad(f.x.options.titleAngle),g=Flotr.getBestTextAli"\
"gn(g.angle,g),Flotr.drawText(e,f.x.options.title,this.plotOffset.left+this.plotWidth/2,this.plotOffset.to"\
"p+f.x.maxLabel.height+this.plotHeight+2*d,g)),f.x2.options.title&&f.x2.used&&(g.textAlign=f.x2.options.ti"\
"tleAlign||\"center\",g.textBaseline=\"bottom\",g.angle=Flotr.toRad(f.x2.options.titleAngle),g=Flotr.getBe"\
"stTextAlign(g.angle,g),Flotr.drawText(e,f.x2.options.title,this.plotOffset.left+this.plotWidth/2,this.plo"\
"tOffset.top-f.x2.maxLabel.height-2*d,g)),f.y.options.title&&f.y.used&&(g.textAlign=f.y.options.titleAlign"\
"||\"right\",g.textBaseline=\"middle\",g.angle=Flotr.toRad(f.y.options.titleAngle),g=Flotr.getBestTextAlig"\
"n(g.angle,g),Flotr.drawText(e,f.y.options.title,this.plotOffset.left-f.y.maxLabel.width-2*d,this.plotOffs"\
"et.top+this.plotHeight/2,g)),f.y2.options.title&&f.y2.used&&(g.textAlign=f.y2.options.titleAlign||\"left\"\
"",g.textBaseline=\"middle\",g.angle=Flotr.toRad(f.y2.options.titleAngle),g=Flotr.getBestTextAlign(g.angle"\
",g),Flotr.drawText(e,f.y2.options.title,this.plotOffset.left+this.plotWidth+f.y2.maxLabel.width+2*d,this."\
"plotOffset.top+this.plotHeight/2,g))}else{b=[],c.title&&b.push('<div style=\"position:absolute;top:0;left"\
":',this.plotOffset.left,\"px;font-size:1em;font-weight:bold;text-align:center;width:\",this.plotWidth,'px"\
";\" class=\"flotr-title\">',c.title,\"</div>\"),c.subtitle&&b.push('<div style=\"position:absolute;top:',"\
"this.titleHeight,\"px;left:\",this.plotOffset.left,\"px;font-size:smaller;text-align:center;width:\",this"\
".plotWidth,'px;\" class=\"flotr-subtitle\">',c.subtitle,\"</div>\"),b.push(\"</div>\"),b.push('<div class"\
"=\"flotr-axis-title\" style=\"font-weight:bold;\">'),f.x.options.title&&f.x.used&&b.push('<div style=\"po"\
"sition:absolute;top:',this.plotOffset.top+this.plotHeight+c.grid.labelMargin+f.x.titleSize.height,\"px;le"\
"ft:\",this.plotOffset.left,\"px;width:\",this.plotWidth,\"px;text-align:\",f.x.options.titleAlign,';\" cl"\
"ass=\"flotr-axis-title flotr-axis-title-x1\">',f.x.options.title,\"</div>\"),f.x2.options.title&&f.x2.use"\
"d&&b.push('<div style=\"position:absolute;top:0;left:',this.plotOffset.left,\"px;width:\",this.plotWidth,"\
"\"px;text-align:\",f.x2.options.titleAlign,';\" class=\"flotr-axis-title flotr-axis-title-x2\">',f.x2.opt"\
"ions.title,\"</div>\"),f.y.options.title&&f.y.used&&b.push('<div style=\"position:absolute;top:',this.plo"\
"tOffset.top+this.plotHeight/2-f.y.titleSize.height/2,\"px;left:0;text-align:\",f.y.options.titleAlign,';\"\
"" class=\"flotr-axis-title flotr-axis-title-y1\">',f.y.options.title,\"</div>\"),f.y2.options.title&&f.y2"\
".used&&b.push('<div style=\"position:absolute;top:',this.plotOffset.top+this.plotHeight/2-f.y.titleSize.h"\
"eight/2,\"px;right:0;text-align:\",f.y2.options.titleAlign,';\" class=\"flotr-axis-title flotr-axis-title"\
"-y2\">',f.y2.options.title,\"</div>\"),b=b.join(\"\");var h=a.create(\"div\");a.setStyles({color:c.grid.c"\
"olor}),h.className=\"flotr-titles\",a.insert(this.el,h),a.insert(h,b)}}})}();\n"\
;

	fprintf(output,"%s\n",flotr2);

}

void barplot(char* figure, char* data, FILE* output){

	fprintf(output,"<div id=\"%s\" style=\"width:800px; height:584px; margin: 8px auto;\"></div>\n", figure);
	fprintf(output,"<script type=\"text/javascript\">\n");
	fprintf(output,"(function basic_bars (container, horizontal) {\n");

	fprintf(output,"  var\n");
	fprintf(output,"    horizontal = (horizontal ? true : false), // Show horizontal bars\n");
	fprintf(output,"    data = [%s];\n", data);
              
  // Draw the graph
	fprintf(output,"  graph = Flotr.draw(container, [ data ], { \n");
	//fprintf(output,"     title: 'Download Image Example',\n");
	fprintf(output,"     bars : { show : true, horizontal : horizontal,\n");
	fprintf(output,"            shadowSize : 0, barWidth : 0.9 },\n");
	fprintf(output,"     mouse: { track : true, relative : true },\n");
	fprintf(output,"     xaxis   : { noTicks : 20 },\n");
	fprintf(output,"     yaxis: { min : 0, autoscaleMargin : 0 }\n");
	fprintf(output,"  });\n");
	fprintf(output,"})(document.getElementById(\"%s\"));\n", figure);
	fprintf(output,"</script>\n");

}

void boxplot(char* figure, char* data, FILE* output){

	fprintf(output,"<div id=\"%s\" style=\"width:800px; height:584px; margin: 8px auto;\"></div>\n", figure);
	fprintf(output,"<script type=\"text/javascript\">\n");
	fprintf(output,"(function basic_candle (container) {\n");

	fprintf(output,"  var d1 = [%s];\n",data);
    
	fprintf(output,"  graph = Flotr.draw(container, [ d1 ], { \n");
	fprintf(output,"     candles : { show : true, candleWidth : 0.6, color : '#000'},\n");
	fprintf(output,"     xaxis   : { noTicks : 10 },\n");
	fprintf(output,"     yaxis : { max : 50, min: 0 }\n");
	fprintf(output,"  });\n");
	fprintf(output,"})(document.getElementById(\"%s\"));\n", figure);
	fprintf(output,"</script>\n");

}

void head(FILE* output){
	fprintf(output,"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n");
	fprintf(output,"<html xmlns=\"http://www.w3.org/1999/xhtml\">\n");
	fprintf(output,"<head>\n");
	fprintf(output,"<script type=\"text/javascript\">");
	flotr2_js(output);
	fprintf(output,"</script>\n");
	fprintf(output,"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n");
	fprintf(output,"<title>Statistics</title>\n");

	fprintf(output,"<style type=\"text/css\">\n");
	fprintf(output,"	body {\n");
	fprintf(output,"        margin: 0px;\n");
	fprintf(output,"        padding: 0px;\n");
	fprintf(output,"	}\n");
	fprintf(output,"	#container {\n");
	fprintf(output,"        width : 800px;\n");
	fprintf(output,"        height: 584px;\n");
	fprintf(output,"        margin: 8px auto;\n");
	fprintf(output,"	}\n");
	fprintf(output,"	table {\n");
	fprintf(output,"	  border-color: #000;\n");
	fprintf(output,"	  border-spacing: 1px;\n");
	fprintf(output,"	  border-style: solid;\n");
	fprintf(output,"	  border-width: 1px;\n");
	fprintf(output,"	  cell-spacing: 2px;\n");
	fprintf(output,"	  background-color: #fff;\n");
	fprintf(output,"	}\n");
	fprintf(output,"	td, th {\n");
	fprintf(output,"	font-family: Arial, Helvetica, sans-serif;\n");
	fprintf(output,"	font-size: 10pt;\n");
	fprintf(output,"	padding: 5px 0.5em;\n");
	fprintf(output,"	white-space: nowrap;\n");
	fprintf(output,"	}\n");
	fprintf(output,"    </style>\n");
	
	fprintf(output,"</head>\n<body>\n");
}
void foot(FILE* output){
	fprintf(output,"</body>\n</html>\n");
}

void plot_cycle(FILE* output, cycle *cycles){

	int i;
	int Q1,Q3,IQR;
	int LeftWisker, RightWisker;
	int total;

	char vectorA[2000];
	char vectorC[2000];
	char vectorG[2000];
	char vectorT[2000];
	char vectorN[2000];

	fprintf(output,"<div>\n<hr>\n<center><b>Statistic of Read Cycles</b></center><p></div>\n");
	fprintf(output,"<div id=\"QBoxplot\" style=\"width:800px; height:584px; margin: 8px auto;\"></div>\n");
	fprintf(output,"<script type=\"text/javascript\">\n");
	fprintf(output,"(function basic_candle (container) {\n");

	fprintf(output,"  var d1 = [");

	for (i=0;i<MAX_READ_LENGTH;i++) {
		if(cycles[i].count == 0)
			continue;
	
		total = cycles[i].nucleotide_count[A] + cycles[i].nucleotide_count[C]
                         + cycles[i].nucleotide_count[G] + cycles[i].nucleotide_count[T]
                         + cycles[i].nucleotide_count[N];
		Q1 = quartile_value ( cycles[i].count / 4, cycles[i].bases_quality_count );
		Q3 = quartile_value ( cycles[i].count * 3 / 4, cycles[i].bases_quality_count );
		IQR = Q3 - Q1 ;
	
		if ( (Q1 - IQR*3/2) < cycles[i].min )
			LeftWisker = cycles[i].min;
		else
			LeftWisker = (Q1 - IQR*3/2); //TODO - make sure there's an observed value at this point
	
		if ( (Q3 + IQR*3/2) > cycles[i].max )
			RightWisker = cycles[i].max;
		else
			RightWisker = (Q3 + IQR*3/2); //TODO - make sure there's an observed value at this point

		if(i == 0){
			fprintf(output,"[%d,%d,%d,%d,%d]", i+1, Q3, RightWisker, LeftWisker, Q1);
			sprintf(vectorA, "[%d,%d]", i+1, cycles[i].nucleotide_count[A]);
			sprintf(vectorC, "[%d,%d]", i+1, cycles[i].nucleotide_count[C]);
			sprintf(vectorG, "[%d,%d]", i+1, cycles[i].nucleotide_count[G]);
			sprintf(vectorT, "[%d,%d]", i+1, cycles[i].nucleotide_count[T]);
			sprintf(vectorN, "[%d,%d]", i+1, cycles[i].nucleotide_count[N]);
		}else{
			fprintf(output,",[%d,%d,%d,%d,%d]", i+1, Q3, RightWisker, LeftWisker, Q1);
			sprintf(vectorA, "%s,[%d,%d]", vectorA, i+1, cycles[i].nucleotide_count[A]);
			sprintf(vectorC, "%s,[%d,%d]", vectorC, i+1, cycles[i].nucleotide_count[C]);
			sprintf(vectorG, "%s,[%d,%d]", vectorG, i+1, cycles[i].nucleotide_count[G]);
			sprintf(vectorT, "%s,[%d,%d]", vectorT, i+1, cycles[i].nucleotide_count[T]);
			sprintf(vectorN, "%s,[%d,%d]", vectorN, i+1, cycles[i].nucleotide_count[N]);
		}
			
	}

	fprintf(output,"];\n"); 
	fprintf(output,"  graph = Flotr.draw(container, [ d1 ], { \n");
	fprintf(output,"     candles : { show : true, candleWidth : 0.6, color : '#000'},\n");
	fprintf(output,"     xaxis   : { noTicks : 20 },\n");
	fprintf(output,"     yaxis : { max : 50, min: 0 }\n");
	fprintf(output,"  });\n");
	fprintf(output,"})(document.getElementById(\"QBoxplot\"));\n");
	fprintf(output,"</script>\n");

	fprintf(output,"<div>\n<hr>\n<center><b>Nucleotide Statistic of Read Cycles</b></center><p></div>\n");
	fprintf(output,"<div id=\"NBarplot\" style=\"width:800px; height:584px; margin: 8px auto;\"></div>\n");
	fprintf(output,"<script type=\"text/javascript\">\n");
	fprintf(output,"(function bars_stacked (container, horizontal) {\n");

	fprintf(output,"  var d1 = [%s],\n", vectorA);
	fprintf(output,"      d2 = [%s],\n", vectorC);
	fprintf(output,"      d3 = [%s],\n", vectorG);
	fprintf(output,"      d4 = [%s],\n", vectorT);
	fprintf(output,"      d5 = [%s];\n", vectorN);

	fprintf(output,"  graph = Flotr.draw(container, [ \n");
	fprintf(output,"    { data : d1, label : 'A' },\n");
	fprintf(output,"    { data : d2, label : 'C' },\n");
	fprintf(output,"    { data : d3, label : 'G' },\n");
	fprintf(output,"    { data : d4, label : 'T' },\n");
	fprintf(output,"    { data : d5, label : 'N' }],\n");
	fprintf(output,"    { \n");
	fprintf(output,"    legend : { noColumns: 5, backgroundColor : '#FFFFFF' },\n");
	fprintf(output,"    bars : { show : true, stacked : true,\n");
	fprintf(output,"      horizontal : horizontal, barWidth : 0.6,\n");
	fprintf(output,"      lineWidth : 1, shadowSize : 0 },\n");
	fprintf(output,"     xaxis   : { noTicks : 20 },\n");
	fprintf(output,"     yaxis   : { max : %d },\n", total);
	fprintf(output,"    grid : { verticalLines : horizontal,\n");
	fprintf(output,"      horizontalLines : !horizontal }\n");
	fprintf(output,"  });\n");
	fprintf(output,"})(document.getElementById(\"NBarplot\"));\n");
	fprintf(output,"</script>\n");

}

void print_stats(FILE* output, stats *seq_stats){

	fprintf(output,"<div>\n<hr>\n<center><b>Statistic Summary</b></center><p>\n");
	fprintf(output,"<table border=\"0\" align=\"center\" width=\"750\">\n");
	fprintf(output,"<tr bgcolor=\"#dddddd\">\n"); 
	fprintf(output,"	<td valign=\"top\" widt=\"280\"> <b> Total read number </b> </td>\n");
	fprintf(output,"	<td width=\"470\"> %ld </td>\n	</tr>\n", seq_stats->total_reads);
	fprintf(output,"<tr bgcolor=\"#ffffff\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Total base number </b> </td>\n");
	fprintf(output,"	<td> %ld </td>\n	</tr>\n", seq_stats->total_length);
	fprintf(output,"<tr bgcolor=\"#dddddd\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Minimum read length </b> </td>\n");
	fprintf(output,"	<td> %d </td>\n	</tr>\n", seq_stats->min_length);
	fprintf(output,"<tr bgcolor=\"#ffffff\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Maximum read length </b> </td>\n");
	fprintf(output,"	<td> %d </td>\n	</tr>\n", seq_stats->max_length);
	fprintf(output,"<tr bgcolor=\"#dddddd\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Mean read length </b> </td>\n");
	fprintf(output,"	<td> %2.2f </td>\n	</tr>\n", seq_stats->mean_length);
	fprintf(output,"<tr bgcolor=\"#ffffff\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Minimum base quality </b> </td>\n");
	fprintf(output,"	<td> %d </td>\n	</tr>\n", seq_stats->min_qual);
	fprintf(output,"<tr bgcolor=\"#dddddd\">\n"); 
	fprintf(output,"	<td valign=\"top\"> <b> Maximum base quality </b> </td>\n");
	fprintf(output,"	<td> %d </td>\n	</tr>\n", seq_stats->max_qual);
	fprintf(output,"</table>\n</div>\n");

	fprintf(output,"<div>\n<hr>\n<center><b>Statistic of Single Base Quality</b></center><p></div>\n");

	int i = 0;
	char vector[BUFFER_LENGTH] = "";
	for(i=0; i<50; i++){
		if(i < seq_stats->min_qual || i > seq_stats->max_qual)
			continue;
		if(strlen(vector) == 0)
			sprintf(vector, "[%d,%ld]", i, seq_stats->qual_hist[i]);
		else
			sprintf(vector, "%s,[%d,%ld]", vector, i, seq_stats->qual_hist[i]);
	}
	barplot("SQbarplot", vector, output);

	fprintf(output,"<div>\n<hr>\n<center><b>Statistic of GC Content of Reads</b></center><p></div>\n");
	for(i=0; i<20; i++){
		if(i == 0)
			sprintf(vector, "[%d,%ld]", i, seq_stats->gc_content[i]);
		else
			sprintf(vector, "%s,[%d,%ld]", vector, i, seq_stats->gc_content[i]);
	}

	barplot("GCbarplot", vector, output);

	if(seq_stats->max_length != seq_stats->min_length){
               	warning_msg("read length is not same, discard cycles statistic information\n");
	}else{
		//boxplot("CyclesBoxplot", vectorQ);
		plot_cycle(output, seq_stats->cycles);
		//stack_barplot("SBarplot");
	}
}
