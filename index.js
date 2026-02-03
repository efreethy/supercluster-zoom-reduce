// Supercluster with zoom-aware reduce

import KDBush from 'kdbush';

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

export default class Supercluster {
    constructor(options) {
        this.options = Object.assign(Object.create(defaultOptions), options || {});

        // Wrap user's reduce to add internal zoom-aware tracking
        const userReduce = this.options.reduce;
        if (userReduce) {
            this._userReduce = userReduce;
            this.options.reduce = (acc, props, zoom) => {
                // Call user's reduce first (may or may not use zoom param)
                userReduce(acc, props, zoom);

                // Internal zoom-aware tracking
                acc._byZoom = acc._byZoom || Object.create(null);
                acc._originZoom = acc._originZoom != null ? acc._originZoom : zoom;

                // Convert propertyStatuses (Set or Array) to array for this zoom bucket
                const statuses = acc.propertyStatuses;
                const statusArray = statuses instanceof Set
                    ? Array.from(statuses)
                    : (Array.isArray(statuses) ? statuses : []);

                // Update zoom bucket with current statuses
                acc._byZoom[zoom] = statusArray.slice(); // copy

                return acc;
            };
        }

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

                            // Initialize zoom tracking and seed with origin point's statuses
                            clusterProperties._byZoom = Object.create(null);
                            clusterProperties._originZoom = zoom;

                            // Seed bucket with origin point's statuses (handles Set or Array)
                            const originStatuses = clusterProperties.propertyStatuses;
                            clusterProperties._byZoom[zoom] = originStatuses instanceof Set
                                ? Array.from(originStatuses)
                                : (Array.isArray(originStatuses) ? originStatuses.slice() : []);

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
    const count = data[i + OFFSET_NUM];
    const abbrev =
        count >= 10000 ? (Math.round(count / 1000) + 'k') :
        count >= 1000 ? ((Math.round(count / 100) / 10) + 'k') : count;
    const propIndex = data[i + OFFSET_PROP];
    const properties = propIndex === -1 ? {} : Object.assign({}, clusterProps[propIndex]);

    // Get statuses from zoom bucket, converting Set to Array for serialization
    let statuses;
    if (properties._byZoom && properties._originZoom != null) {
        statuses = properties._byZoom[properties._originZoom];
    } else if (properties.propertyStatuses) {
        statuses = properties.propertyStatuses instanceof Set
            ? Array.from(properties.propertyStatuses)
            : properties.propertyStatuses;
    }

    const extra = {
        cluster: true,
        cluster_id: data[i + OFFSET_ID],
        point_count: count,
        point_count_abbreviated: abbrev
    };

    if (statuses && statuses.length > 0) {
        extra.propertyStatuses = statuses;
    }

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
