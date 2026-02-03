(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
typeof define === 'function' && define.amd ? define(factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, global.Supercluster = factory());
})(this, (function () { 'use strict';

const ARRAY_TYPES = [
    Int8Array, Uint8Array, Uint8ClampedArray, Int16Array, Uint16Array,
    Int32Array, Uint32Array, Float32Array, Float64Array
];

/** @typedef {Int8ArrayConstructor | Uint8ArrayConstructor | Uint8ClampedArrayConstructor | Int16ArrayConstructor | Uint16ArrayConstructor | Int32ArrayConstructor | Uint32ArrayConstructor | Float32ArrayConstructor | Float64ArrayConstructor} TypedArrayConstructor */

const VERSION = 1; // serialized format version
const HEADER_SIZE = 8;

class KDBush {

    /**
     * Creates an index from raw `ArrayBuffer` data.
     * @param {ArrayBuffer} data
     */
    static from(data) {
        if (!(data instanceof ArrayBuffer)) {
            throw new Error('Data must be an instance of ArrayBuffer.');
        }
        const [magic, versionAndType] = new Uint8Array(data, 0, 2);
        if (magic !== 0xdb) {
            throw new Error('Data does not appear to be in a KDBush format.');
        }
        const version = versionAndType >> 4;
        if (version !== VERSION) {
            throw new Error(`Got v${version} data when expected v${VERSION}.`);
        }
        const ArrayType = ARRAY_TYPES[versionAndType & 0x0f];
        if (!ArrayType) {
            throw new Error('Unrecognized array type.');
        }
        const [nodeSize] = new Uint16Array(data, 2, 1);
        const [numItems] = new Uint32Array(data, 4, 1);

        return new KDBush(numItems, nodeSize, ArrayType, data);
    }

    /**
     * Creates an index that will hold a given number of items.
     * @param {number} numItems
     * @param {number} [nodeSize=64] Size of the KD-tree node (64 by default).
     * @param {TypedArrayConstructor} [ArrayType=Float64Array] The array type used for coordinates storage (`Float64Array` by default).
     * @param {ArrayBuffer} [data] (For internal use only)
     */
    constructor(numItems, nodeSize = 64, ArrayType = Float64Array, data) {
        if (isNaN(numItems) || numItems < 0) throw new Error(`Unpexpected numItems value: ${numItems}.`);

        this.numItems = +numItems;
        this.nodeSize = Math.min(Math.max(+nodeSize, 2), 65535);
        this.ArrayType = ArrayType;
        this.IndexArrayType = numItems < 65536 ? Uint16Array : Uint32Array;

        const arrayTypeIndex = ARRAY_TYPES.indexOf(this.ArrayType);
        const coordsByteSize = numItems * 2 * this.ArrayType.BYTES_PER_ELEMENT;
        const idsByteSize = numItems * this.IndexArrayType.BYTES_PER_ELEMENT;
        const padCoords = (8 - idsByteSize % 8) % 8;

        if (arrayTypeIndex < 0) {
            throw new Error(`Unexpected typed array class: ${ArrayType}.`);
        }

        if (data && (data instanceof ArrayBuffer)) { // reconstruct an index from a buffer
            this.data = data;
            this.ids = new this.IndexArrayType(this.data, HEADER_SIZE, numItems);
            this.coords = new this.ArrayType(this.data, HEADER_SIZE + idsByteSize + padCoords, numItems * 2);
            this._pos = numItems * 2;
            this._finished = true;
        } else { // initialize a new index
            this.data = new ArrayBuffer(HEADER_SIZE + coordsByteSize + idsByteSize + padCoords);
            this.ids = new this.IndexArrayType(this.data, HEADER_SIZE, numItems);
            this.coords = new this.ArrayType(this.data, HEADER_SIZE + idsByteSize + padCoords, numItems * 2);
            this._pos = 0;
            this._finished = false;

            // set header
            new Uint8Array(this.data, 0, 2).set([0xdb, (VERSION << 4) + arrayTypeIndex]);
            new Uint16Array(this.data, 2, 1)[0] = nodeSize;
            new Uint32Array(this.data, 4, 1)[0] = numItems;
        }
    }

    /**
     * Add a point to the index.
     * @param {number} x
     * @param {number} y
     * @returns {number} An incremental index associated with the added item (starting from `0`).
     */
    add(x, y) {
        const index = this._pos >> 1;
        this.ids[index] = index;
        this.coords[this._pos++] = x;
        this.coords[this._pos++] = y;
        return index;
    }

    /**
     * Perform indexing of the added points.
     */
    finish() {
        const numAdded = this._pos >> 1;
        if (numAdded !== this.numItems) {
            throw new Error(`Added ${numAdded} items when expected ${this.numItems}.`);
        }
        // kd-sort both arrays for efficient search
        sort(this.ids, this.coords, this.nodeSize, 0, this.numItems - 1, 0);

        this._finished = true;
        return this;
    }

    /**
     * Search the index for items within a given bounding box.
     * @param {number} minX
     * @param {number} minY
     * @param {number} maxX
     * @param {number} maxY
     * @returns {number[]} An array of indices correponding to the found items.
     */
    range(minX, minY, maxX, maxY) {
        if (!this._finished) throw new Error('Data not yet indexed - call index.finish().');

        const {ids, coords, nodeSize} = this;
        const stack = [0, ids.length - 1, 0];
        const result = [];

        // recursively search for items in range in the kd-sorted arrays
        while (stack.length) {
            const axis = stack.pop() || 0;
            const right = stack.pop() || 0;
            const left = stack.pop() || 0;

            // if we reached "tree node", search linearly
            if (right - left <= nodeSize) {
                for (let i = left; i <= right; i++) {
                    const x = coords[2 * i];
                    const y = coords[2 * i + 1];
                    if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push(ids[i]);
                }
                continue;
            }

            // otherwise find the middle index
            const m = (left + right) >> 1;

            // include the middle item if it's in range
            const x = coords[2 * m];
            const y = coords[2 * m + 1];
            if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push(ids[m]);

            // queue search in halves that intersect the query
            if (axis === 0 ? minX <= x : minY <= y) {
                stack.push(left);
                stack.push(m - 1);
                stack.push(1 - axis);
            }
            if (axis === 0 ? maxX >= x : maxY >= y) {
                stack.push(m + 1);
                stack.push(right);
                stack.push(1 - axis);
            }
        }

        return result;
    }

    /**
     * Search the index for items within a given radius.
     * @param {number} qx
     * @param {number} qy
     * @param {number} r Query radius.
     * @returns {number[]} An array of indices correponding to the found items.
     */
    within(qx, qy, r) {
        if (!this._finished) throw new Error('Data not yet indexed - call index.finish().');

        const {ids, coords, nodeSize} = this;
        const stack = [0, ids.length - 1, 0];
        const result = [];
        const r2 = r * r;

        // recursively search for items within radius in the kd-sorted arrays
        while (stack.length) {
            const axis = stack.pop() || 0;
            const right = stack.pop() || 0;
            const left = stack.pop() || 0;

            // if we reached "tree node", search linearly
            if (right - left <= nodeSize) {
                for (let i = left; i <= right; i++) {
                    if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2) result.push(ids[i]);
                }
                continue;
            }

            // otherwise find the middle index
            const m = (left + right) >> 1;

            // include the middle item if it's in range
            const x = coords[2 * m];
            const y = coords[2 * m + 1];
            if (sqDist(x, y, qx, qy) <= r2) result.push(ids[m]);

            // queue search in halves that intersect the query
            if (axis === 0 ? qx - r <= x : qy - r <= y) {
                stack.push(left);
                stack.push(m - 1);
                stack.push(1 - axis);
            }
            if (axis === 0 ? qx + r >= x : qy + r >= y) {
                stack.push(m + 1);
                stack.push(right);
                stack.push(1 - axis);
            }
        }

        return result;
    }
}

/**
 * @param {Uint16Array | Uint32Array} ids
 * @param {InstanceType<TypedArrayConstructor>} coords
 * @param {number} nodeSize
 * @param {number} left
 * @param {number} right
 * @param {number} axis
 */
function sort(ids, coords, nodeSize, left, right, axis) {
    if (right - left <= nodeSize) return;

    const m = (left + right) >> 1; // middle index

    // sort ids and coords around the middle index so that the halves lie
    // either left/right or top/bottom correspondingly (taking turns)
    select(ids, coords, m, left, right, axis);

    // recursively kd-sort first half and second half on the opposite axis
    sort(ids, coords, nodeSize, left, m - 1, 1 - axis);
    sort(ids, coords, nodeSize, m + 1, right, 1 - axis);
}

/**
 * Custom Floyd-Rivest selection algorithm: sort ids and coords so that
 * [left..k-1] items are smaller than k-th item (on either x or y axis)
 * @param {Uint16Array | Uint32Array} ids
 * @param {InstanceType<TypedArrayConstructor>} coords
 * @param {number} k
 * @param {number} left
 * @param {number} right
 * @param {number} axis
 */
function select(ids, coords, k, left, right, axis) {

    while (right > left) {
        if (right - left > 600) {
            const n = right - left + 1;
            const m = k - left + 1;
            const z = Math.log(n);
            const s = 0.5 * Math.exp(2 * z / 3);
            const sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            const newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            const newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            select(ids, coords, k, newLeft, newRight, axis);
        }

        const t = coords[2 * k + axis];
        let i = left;
        let j = right;

        swapItem(ids, coords, left, k);
        if (coords[2 * right + axis] > t) swapItem(ids, coords, left, right);

        while (i < j) {
            swapItem(ids, coords, i, j);
            i++;
            j--;
            while (coords[2 * i + axis] < t) i++;
            while (coords[2 * j + axis] > t) j--;
        }

        if (coords[2 * left + axis] === t) swapItem(ids, coords, left, j);
        else {
            j++;
            swapItem(ids, coords, j, right);
        }

        if (j <= k) left = j + 1;
        if (k <= j) right = j - 1;
    }
}

/**
 * @param {Uint16Array | Uint32Array} ids
 * @param {InstanceType<TypedArrayConstructor>} coords
 * @param {number} i
 * @param {number} j
 */
function swapItem(ids, coords, i, j) {
    swap(ids, i, j);
    swap(coords, 2 * i, 2 * j);
    swap(coords, 2 * i + 1, 2 * j + 1);
}

/**
 * @param {InstanceType<TypedArrayConstructor>} arr
 * @param {number} i
 * @param {number} j
 */
function swap(arr, i, j) {
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

/**
 * @param {number} ax
 * @param {number} ay
 * @param {number} bx
 * @param {number} by
 */
function sqDist(ax, ay, bx, by) {
    const dx = ax - bx;
    const dy = ay - by;
    return dx * dx + dy * dy;
}

// Supercluster with zoom-aware reduce


const defaultOptions = {
    minZoom: 0,
    maxZoom: 16,
    minPoints: 2,
    radius: 40,
    extent: 512,
    nodeSize: 64,
    log: false,
    generateId: false,
    reduce: null,
    map: function (props) { return props; }
};

const fround = Math.fround || (function (tmp) { return function (x) { tmp[0] = +x; return tmp[0]; }; })(new Float32Array(1));

const OFFSET_ZOOM = 2;
const OFFSET_ID = 3;
const OFFSET_PARENT = 4;
const OFFSET_NUM = 5;
const OFFSET_PROP = 6;

class Supercluster {
    constructor(options) {
        this.options = Object.assign(Object.create(defaultOptions), options || {});
        this.trees = new Array(this.options.maxZoom + 1);
        this.stride = this.options.reduce ? 7 : 6;
        this.clusterProps = [];
        this.points = [];
    }

    load(points) {
        const log = this.options.log;
        const minZoom = this.options.minZoom;
        const maxZoom = this.options.maxZoom;

        if (log) console.time('total time');

        const timerId = `prepare ${points.length} points`;
        if (log) console.time(timerId);

        this.points = points;
        const data = [];

        for (let i = 0; i < points.length; i++) {
            const p = points[i];
            if (!p || !p.geometry) continue;

            const coords = p.geometry.coordinates;
            if (!Array.isArray(coords) || coords.length < 2) continue;
            const lng = coords[0];
            const lat = coords[1];

            const x = fround(lngX(lng));
            const y = fround(latY(lat));
            data.push(x, y, Infinity, i, -1, 1);
            if (this.options.reduce) data.push(0);
        }
        let tree = (this.trees[maxZoom + 1] = this._createTree(data));

        if (log) console.timeEnd(timerId);

        for (let z = maxZoom; z >= minZoom; z--) {
            const now = +Date.now();
            tree = this.trees[z] = this._createTree(this._cluster(tree, z));
            if (log) console.log('z%d: %d clusters in %dms', z, tree.numItems, +Date.now() - now);
        }

        if (log) console.timeEnd('total time');
        return this;
    }

    getClusters(bbox, zoom) {
        let minLng = ((bbox[0] + 180) % 360 + 360) % 360 - 180;
        const minLat = Math.max(-90, Math.min(90, bbox[1]));
        let maxLng = bbox[2] === 180 ? 180 : ((bbox[2] + 180) % 360 + 360) % 360 - 180;
        const maxLat = Math.max(-90, Math.min(90, bbox[3]));

        if (bbox[2] - bbox[0] >= 360) {
            minLng = -180;
            maxLng = 180;
        } else if (minLng > maxLng) {
            const easternHem = this.getClusters([minLng, minLat, 180, maxLat], zoom);
            const westernHem = this.getClusters([-180, minLat, maxLng, maxLat], zoom);
            return easternHem.concat(westernHem);
        }

        const tree = this.trees[this._limitZoom(zoom)];
        const ids = tree.range(lngX(minLng), latY(maxLat), lngX(maxLng), latY(minLat));
        const data = tree.data;
        const clusters = [];
        for (let t = 0; t < ids.length; t++) {
            const id = ids[t];
            const k = this.stride * id;
            clusters.push(
                data[k + OFFSET_NUM] > 1
                    ? getClusterJSON(data, k, this.clusterProps)
                    : this.points[data[k + OFFSET_ID]]
            );
        }
        return clusters;
    }

    getChildren(clusterId) {
        const originId = this._getOriginId(clusterId);
        const originZoom = this._getOriginZoom(clusterId);
        const errorMsg = 'No cluster with the specified id.';

        const tree = this.trees[originZoom];
        if (!tree) throw new Error(errorMsg);

        const data = tree.data;
        if (originId * this.stride >= data.length) throw new Error(errorMsg);

        const r = this.options.radius / (this.options.extent * Math.pow(2, originZoom - 1));
        const x = data[originId * this.stride];
        const y = data[originId * this.stride + 1];
        const ids = tree.within(x, y, r);
        const children = [];
        for (let t = 0; t < ids.length; t++) {
            const id = ids[t];
            const k = id * this.stride;
            if (data[k + OFFSET_PARENT] === clusterId) {
                children.push(
                    data[k + OFFSET_NUM] > 1
                        ? getClusterJSON(data, k, this.clusterProps)
                        : this.points[data[k + OFFSET_ID]]
                );
            }
        }

        if (children.length === 0) throw new Error(errorMsg);
        return children;
    }

    getLeaves(clusterId, limit, offset) {
        limit = limit || 10;
        offset = offset || 0;
        const leaves = [];
        this._appendLeaves(leaves, clusterId, limit, offset, 0);
        return leaves;
    }

    getTile(z, x, y) {
        const tree = this.trees[this._limitZoom(z)];
        const z2 = Math.pow(2, z);
        const extent = this.options.extent;
        const radius = this.options.radius;
        const p = radius / extent;
        const top = (y - p) / z2;
        const bottom = (y + 1 + p) / z2;

        const tile = { features: [] };

        this._addTileFeatures(
            tree.range((x - p) / z2, top, (x + 1 + p) / z2, bottom),
            tree.data, x, y, z2, tile);

        if (x === 0) {
            this._addTileFeatures(
                tree.range(1 - p / z2, top, 1, bottom),
                tree.data, z2, y, z2, tile);
        }
        if (x === z2 - 1) {
            this._addTileFeatures(
                tree.range(0, top, p / z2, bottom),
                tree.data, -1, y, z2, tile);
        }

        return tile.features.length ? tile : null;
    }

    getClusterExpansionZoom(clusterId) {
        let expansionZoom = this._getOriginZoom(clusterId) - 1;
        while (expansionZoom <= this.options.maxZoom) {
            const children = this.getChildren(clusterId);
            expansionZoom++;
            if (children.length !== 1) break;
            clusterId = children[0].properties.cluster_id;
        }
        return expansionZoom;
    }

    _appendLeaves(result, clusterId, limit, offset, skipped) {
        const children = this.getChildren(clusterId);

        for (let c = 0; c < children.length; c++) {
            const child = children[c];
            const props = child.properties;

            if (props && props.cluster) {
                if (skipped + props.point_count <= offset) {
                    skipped += props.point_count;
                } else {
                    skipped = this._appendLeaves(result, props.cluster_id, limit, offset, skipped);
                }
            } else if (skipped < offset) {
                skipped++;
            } else {
                result.push(child);
            }
            if (result.length === limit) break;
        }

        return skipped;
    }

    _createTree(data) {
        const tree = new KDBush((data.length / this.stride) | 0, this.options.nodeSize, Float32Array);
        for (let i = 0; i < data.length; i += this.stride) tree.add(data[i], data[i + 1]);
        tree.finish();
        tree.data = data;
        return tree;
    }

    _addTileFeatures(ids, data, x, y, z2, tile) {
        for (let idx = 0; idx < ids.length; idx++) {
            const i = ids[idx];
            const k = i * this.stride;
            const isCluster = data[k + OFFSET_NUM] > 1;

            let tags, px, py;
            if (isCluster) {
                tags = getClusterProperties(data, k, this.clusterProps);
                px = data[k];
                py = data[k + 1];
            } else {
                const p = this.points[data[k + OFFSET_ID]];
                tags = p.properties;
                const lng = p.geometry.coordinates[0];
                const lat = p.geometry.coordinates[1];
                px = lngX(lng);
                py = latY(lat);
            }

            const f = {
                type: 1,
                geometry: [[
                    Math.round(this.options.extent * (px * z2 - x)),
                    Math.round(this.options.extent * (py * z2 - y))
                ]],
                tags
            };

            let id;
            if (isCluster || this.options.generateId) {
                id = data[k + OFFSET_ID];
            } else {
                id = this.points[data[k + OFFSET_ID]].id;
            }

            if (id !== undefined) f.id = id;
            tile.features.push(f);
        }
    }

    _limitZoom(z) {
        return Math.max(this.options.minZoom, Math.min(Math.floor(+z), this.options.maxZoom + 1));
    }

    _cluster(tree, zoom) {
        const radius = this.options.radius;
        const extent = this.options.extent;
        const reduce = this.options.reduce;
        const minPoints = this.options.minPoints;
        const r = radius / (extent * Math.pow(2, zoom));
        const data = tree.data;
        const nextData = [];
        const stride = this.stride;

        for (let i = 0; i < data.length; i += stride) {
            if (data[i + OFFSET_ZOOM] <= zoom) continue;
            data[i + OFFSET_ZOOM] = zoom;

            const x = data[i];
            const y = data[i + 1];
            const neighborIds = tree.within(data[i], data[i + 1], r);

            const numPointsOrigin = data[i + OFFSET_NUM];
            let numPoints = numPointsOrigin;

            for (let ni = 0; ni < neighborIds.length; ni++) {
                const neighborId = neighborIds[ni];
                const k = neighborId * stride;
                if (data[k + OFFSET_ZOOM] > zoom) numPoints += data[k + OFFSET_NUM];
            }

            if (numPoints > numPointsOrigin && numPoints >= minPoints) {
                let wx = x * numPointsOrigin;
                let wy = y * numPointsOrigin;

                let clusterProperties;
                let clusterPropIndex = -1;

                const id = (((i / stride) | 0) << 5) + (zoom + 1) + this.points.length;

                for (let ni = 0; ni < neighborIds.length; ni++) {
                    const neighborId = neighborIds[ni];
                    const k = neighborId * stride;

                    if (data[k + OFFSET_ZOOM] <= zoom) continue;
                    data[k + OFFSET_ZOOM] = zoom;

                    const numPoints2 = data[k + OFFSET_NUM];
                    wx += data[k] * numPoints2;
                    wy += data[k + 1] * numPoints2;

                    data[k + OFFSET_PARENT] = id;

                    if (reduce) {
                        if (!clusterProperties) {
                            clusterProperties = this._map(data, i, true);

                            if (!clusterProperties._byZoom) clusterProperties._byZoom = Object.create(null);
                            // Always set _originZoom to current zoom level
                            clusterProperties._originZoom = zoom;
                            if (!clusterProperties._byZoom[zoom]) clusterProperties._byZoom[zoom] = [];

                            clusterPropIndex = this.clusterProps.length;
                            this.clusterProps.push(clusterProperties);
                        }
                        reduce(clusterProperties, this._map(data, k), zoom);
                    }
                }

                data[i + OFFSET_PARENT] = id;
                nextData.push(wx / numPoints, wy / numPoints, Infinity, id, -1, numPoints);
                if (reduce) nextData.push(clusterPropIndex);

            } else {
                for (let j = 0; j < stride; j++) nextData.push(data[i + j]);

                if (numPoints > 1) {
                    for (let ni = 0; ni < neighborIds.length; ni++) {
                        const neighborId = neighborIds[ni];
                        const k = neighborId * stride;
                        if (data[k + OFFSET_ZOOM] <= zoom) continue;
                        data[k + OFFSET_ZOOM] = zoom;
                        for (let j = 0; j < stride; j++) nextData.push(data[k + j]);
                    }
                }
            }
        }

        return nextData;
    }

    _getOriginId(clusterId) {
        return (clusterId - this.points.length) >> 5;
    }

    _getOriginZoom(clusterId) {
        return (clusterId - this.points.length) % 32;
    }

    _map(data, i, clone) {
        if (data[i + OFFSET_NUM] > 1) {
            const props = this.clusterProps[data[i + OFFSET_PROP]];
            return clone ? Object.assign({}, props) : props;
        }
        const original = this.points[data[i + OFFSET_ID]].properties;
        const result = this.options.map(original);
        return clone && result === original ? Object.assign({}, result) : result;
    }
}

function getClusterJSON(data, i, clusterProps) {
    return {
        type: 'Feature',
        id: data[i + OFFSET_ID],
        properties: getClusterProperties(data, i, clusterProps),
        geometry: {
            type: 'Point',
            coordinates: [xLng(data[i]), yLat(data[i + 1])]
        }
    };
}

function getClusterProperties(data, i, clusterProps) {
    var count = data[i + OFFSET_NUM];
    var abbrev =
        count >= 10000 ? (Math.round(count / 1000) + 'k') :
        count >= 1000 ? ((Math.round(count / 100) / 10) + 'k') : count;
    var propIndex = data[i + OFFSET_PROP];
    var properties = propIndex === -1 ? {} : Object.assign({}, clusterProps[propIndex]);

    // Get statuses - try _byZoom[_originZoom] first, then fall back to propertyStatuses
    var byZoom = properties._byZoom || null;
    var originZoom = properties._originZoom;
    var statuses =
        byZoom && originZoom != null && Array.isArray(byZoom[originZoom])
            ? byZoom[originZoom]
            : (Array.isArray(properties.propertyStatuses) ? properties.propertyStatuses : undefined);

    var extra = {
        cluster: true,
        'cluster_id': data[i + OFFSET_ID],
        'point_count': count,
        'point_count_abbreviated': abbrev
    };
    if (statuses && statuses.length > 0) extra.propertyStatuses = statuses;

    return Object.assign(properties, extra);
}

function lngX(lng) {
    return lng / 360 + 0.5;
}

function latY(lat) {
    const sin = Math.sin(lat * Math.PI / 180);
    const y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
    return y < 0 ? 0 : y > 1 ? 1 : y;
}

function xLng(x) {
    return (x - 0.5) * 360;
}

function yLat(y) {
    const y2 = (180 - y * 360) * Math.PI / 180;
    return (360 * Math.atan(Math.exp(y2)) / Math.PI) - 90;
}

return Supercluster;

}));
