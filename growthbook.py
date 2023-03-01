#!/usr/bin/env python
u"""
This is the Python client library for GrowthBook, the open-source
feature flagging and A/B testing platform.
More info at https://www.growthbook.io
"""

from __future__ import division
from __future__ import absolute_import
import re
from urlparse import urlparse, parse_qs


def fnv1a32(unicode):
    hval = 0x811C9DC5
    prime = 0x01000193
    uint32_max = 2 ** 32
    for s in unicode:
        hval = hval ^ ord(s)
        hval = (hval * prime) % uint32_max
    return hval


def gbhash(unicode):
    n = fnv1a32(unicode)
    return (n % 1000) / 1000


def inNamespace(userId, namespace):
    n = gbhash(userId + u"__" + namespace[0])
    return n >= namespace[1] and n < namespace[2]


def getEqualWeights(numVariations):
    if numVariations < 1:
        return []
    return [1 / numVariations for i in xrange(numVariations)]


def getBucketRanges(
    numVariations, coverage = 1, weights = None
):
    if coverage < 0:
        coverage = 0
    if coverage > 1:
        coverage = 1
    if weights is None:
        weights = getEqualWeights(numVariations)
    if len(weights) != numVariations:
        weights = getEqualWeights(numVariations)
    if sum(weights) < 0.99 or sum(weights) > 1.01:
        weights = getEqualWeights(numVariations)

    cumulative = 0
    ranges = []
    for w in weights:
        start = cumulative
        cumulative += w
        ranges.append((start, start + coverage * w))

    return ranges


def chooseVariation(n, ranges):
    for i, r in enumerate(ranges):
        if n >= r[0] and n < r[1]:
            return i
    return -1


def getQueryStringOverride(id, url, numVariations):
    res = urlparse(url)
    if not res.query:
        return None
    qs = parse_qs(res.query)
    if id not in qs:
        return None
    variation = qs[id][0]
    if variation is None or not variation.isdigit():
        return None
    varId = int(variation)
    if varId < 0 or varId >= numVariations:
        return None
    return varId


def evalCondition(attributes, condition):
    if u"$or" in condition:
        return evalOr(attributes, condition[u"$or"])
    if u"$nor" in condition:
        return not evalOr(attributes, condition[u"$nor"])
    if u"$and" in condition:
        return evalAnd(attributes, condition[u"$and"])
    if u"$not" in condition:
        return not evalCondition(attributes, condition[u"$not"])

    for key, value in condition.items():
        if not evalConditionValue(value, getPath(attributes, key)):
            return False

    return True


def evalOr(attributes, conditions):
    if len(conditions) == 0:
        return True

    for condition in conditions:
        if evalCondition(attributes, condition):
            return True
    return False


def evalAnd(attributes, conditions):
    for condition in conditions:
        if not evalCondition(attributes, condition):
            return False
    return True


def isOperatorObject(obj):
    for key in obj.keys():
        if key[0] != u"$":
            return False
    return True


def getType(attributeValue):
    t = type(attributeValue)

    if attributeValue is None:
        return u"null"
    if t is int or t is float:
        return u"number"
    if t is unicode:
        return u"string"
    if t is list or t is set:
        return u"array"
    if t is dict:
        return u"object"
    if t is bool:
        return u"boolean"
    return u"unknown"


def getPath(attributes, path):
    current = attributes
    for segment in path.split(u"."):
        if type(current) is dict and segment in current:
            current = current[segment]
        else:
            return None
    return current


def evalConditionValue(conditionValue, attributeValue):
    if type(conditionValue) is dict and isOperatorObject(conditionValue):
        for key, value in conditionValue.items():
            if not evalOperatorCondition(key, attributeValue, value):
                return False
        return True
    return conditionValue == attributeValue


def elemMatch(condition, attributeValue):
    if not type(attributeValue) is list:
        return False

    for item in attributeValue:
        if isOperatorObject(condition):
            if evalConditionValue(condition, item):
                return True
        else:
            if evalCondition(item, condition):
                return True

    return False


def evalOperatorCondition(operator, attributeValue, conditionValue):
    if operator == u"$eq":
        return attributeValue == conditionValue
    elif operator == u"$ne":
        return attributeValue != conditionValue
    elif operator == u"$lt":
        return attributeValue < conditionValue
    elif operator == u"$lte":
        return attributeValue <= conditionValue
    elif operator == u"$gt":
        return attributeValue > conditionValue
    elif operator == u"$gte":
        return attributeValue >= conditionValue
    elif operator == u"$regex":
        try:
            r = re.compile(conditionValue)
            return bool(r.search(attributeValue))
        except Exception:
            return False
    elif operator == u"$in":
        return attributeValue in conditionValue
    elif operator == u"$nin":
        return not (attributeValue in conditionValue)
    elif operator == u"$elemMatch":
        return elemMatch(conditionValue, attributeValue)
    elif operator == u"$size":
        if not (type(attributeValue) is list):
            return False
        return evalConditionValue(conditionValue, len(attributeValue))
    elif operator == u"$all":
        if not (type(attributeValue) is list):
            return False
        for cond in conditionValue:
            passing = False
            for attr in attributeValue:
                if evalConditionValue(cond, attr):
                    passing = True
            if not passing:
                return False
        return True
    elif operator == u"$exists":
        if not conditionValue:
            return attributeValue is None
        return attributeValue is not None
    elif operator == u"$type":
        return getType(attributeValue) == conditionValue
    elif operator == u"$not":
        return not evalConditionValue(conditionValue, attributeValue)
    return False


class Experiment(object):
    def __init__(
        self,
        key,
        variations,
        weights = None,
        active = True,
        status = u"running",
        coverage = 1,
        condition = None,
        namespace = None,
        url = u"",
        include=None,
        groups = None,
        force = None,
        hashAttribute = u"id",
    ):
        self.key = key
        self.variations = variations
        self.weights = weights
        self.active = active
        self.coverage = coverage
        self.condition = condition
        self.namespace = namespace
        self.force = force
        self.hashAttribute = hashAttribute

        # Deprecated properties
        self.status = status
        self.url = url
        self.include = include
        self.groups = groups

    def to_dict(self):
        return {
            u"key": self.key,
            u"variations": self.variations,
            u"weights": self.weights,
            u"active": self.active,
            u"coverage": self.coverage,
            u"condition": self.condition,
            u"namespace": self.namespace,
            u"force": self.force,
            u"hashAttribute": self.hashAttribute,
        }

    def update(self, data):
        weights = data.get(u"weights", None)
        status = data.get(u"status", None)
        coverage = data.get(u"coverage", None)
        url = data.get(u"url", None)
        groups = data.get(u"groups", None)
        force = data.get(u"force", None)

        if weights is not None:
            self.weights = weights
        if status is not None:
            self.status = status
        if coverage is not None:
            self.coverage = coverage
        if url is not None:
            self.url = url
        if groups is not None:
            self.groups = groups
        if force is not None:
            self.force = force


class Result(object):
    def __init__(
        self,
        variationId,
        inExperiment,
        value,
        hashUsed,
        hashAttribute,
        hashValue,
        featureId,
    ):
        self.variationId = variationId
        self.inExperiment = inExperiment
        self.value = value
        self.hashUsed = hashUsed
        self.hashAttribute = hashAttribute
        self.hashValue = hashValue
        self.featureId = featureId or None

    def to_dict(self):
        return {
            u"featureId": self.featureId,
            u"variationId": self.variationId,
            u"inExperiment": self.inExperiment,
            u"value": self.value,
            u"hashUsed": self.hashUsed,
            u"hashAttribute": self.hashAttribute,
            u"hashValue": self.hashValue,
        }


class Feature(object):
    def __init__(self, defaultValue=None, rules = []):
        self.defaultValue = defaultValue
        self.rules = []
        for rule in rules:
            if isinstance(rule, FeatureRule):
                self.rules.append(rule)
            else:
                self.rules.append(FeatureRule(**rule))

    def to_dict(self):
        return {
            u"defaultValue": self.defaultValue,
            u"rules": [rule.to_dict() for rule in self.rules],
        }


class FeatureRule(object):
    def __init__(
        self,
        key = u"",
        variations = None,
        weights = None,
        coverage = 1,
        condition = None,
        namespace = None,
        force=None,
        hashAttribute = u"id",
    ):
        self.key = key
        self.variations = variations
        self.weights = weights
        self.coverage = coverage
        self.condition = condition
        self.namespace = namespace
        self.force = force
        self.hashAttribute = hashAttribute

    def to_dict(self):
        data = {}
        if self.key:
            data[u"key"] = self.key
        if self.variations is not None:
            data[u"variations"] = self.variations
        if self.weights is not None:
            data[u"weights"] = self.weights
        if self.coverage != 1:
            data[u"coverage"] = self.coverage
        if self.condition is not None:
            data[u"condition"] = self.condition
        if self.namespace is not None:
            data[u"namespace"] = self.namespace
        if self.force is not None:
            data[u"force"] = self.force
        if self.hashAttribute != u"id":
            data[u"hashAttribute"] = self.hashAttribute

        return data


class FeatureResult(object):
    def __init__(
        self,
        value,
        source,
        experiment = None,
        experimentResult = None,
    ):
        self.value = value
        self.source = source
        self.experiment = experiment
        self.experimentResult = experimentResult
        self.on = bool(value)
        self.off = not bool(value)

    def to_dict(self):
        data = {
            u"value": self.value,
            u"source": self.source,
            u"on": self.on,
            u"off": self.off,
        }
        if self.experiment:
            data[u"experiment"] = self.experiment.to_dict()
        if self.experimentResult:
            data[u"experimentResult"] = self.experimentResult.to_dict()

        return data


class GrowthBook(object):
    def __init__(
        self,
        enabled = True,
        attributes = {},
        url = u"",
        features = {},
        qaMode = False,
        trackingCallback=None,
        # Deprecated args
        user = {},
        groups = {},
        overrides = {},
        forcedVariations = {},
    ):
        self._enabled = enabled
        self._attributes = attributes
        self._url = url
        self._features = {}

        if features:
            self.setFeatures(features)

        self._qaMode = qaMode
        self._trackingCallback = trackingCallback

        # Deprecated args
        self._user = user
        self._groups = groups
        self._overrides = overrides
        self._forcedVariations = forcedVariations

        self._tracked = {}
        self._assigned = {}
        self._subscriptions = set()

    def setFeatures(self, features):
        self._features = {}
        for key, feature in features.items():
            if isinstance(feature, Feature):
                self._features[key] = feature
            else:
                self._features[key] = Feature(**feature)

    def getFeatures(self):
        return self._features

    def setAttributes(self, attributes):
        self._attributes = attributes

    def getAttributes(self):
        return self._attributes

    def destroy(self):
        self._subscriptions.clear()
        self._tracked.clear()
        self._assigned.clear()
        self._trackingCallback = None
        self._forcedVariations.clear()
        self._overrides.clear()
        self._groups.clear()
        self._attributes.clear()
        self._features.clear()

    def isOn(self, key):
        return self.evalFeature(key).on

    def isOff(self, key):
        return self.evalFeature(key).off

    def getFeatureValue(self, key, fallback):
        res = self.evalFeature(key)
        return res.value if res.value is not None else fallback

    def evalFeature(self, key):
        if key not in self._features:
            return FeatureResult(None, u"unknownFeature")

        feature = self._features[key]
        for rule in feature.rules:
            if rule.condition:
                if not evalCondition(self._attributes, rule.condition):
                    continue
            if rule.force is not None:
                if rule.coverage < 1:
                    hashValue = self._getHashValue(rule.hashAttribute)
                    if not hashValue:
                        continue

                    n = gbhash(hashValue + key)

                    if n > rule.coverage:
                        continue
                return FeatureResult(rule.force, u"force")

            if rule.variations is None:
                continue

            exp = Experiment(
                key=rule.key or key,
                variations=rule.variations,
                coverage=rule.coverage,
                weights=rule.weights,
                hashAttribute=rule.hashAttribute,
                namespace=rule.namespace,
            )

            result = self._run(exp, key)
            self._fireSubscriptions(exp, result)

            if not result.inExperiment:
                continue

            return FeatureResult(result.value, u"experiment", exp, result)

        return FeatureResult(feature.defaultValue, u"defaultValue")

    def getAllResults(self):
        return self._assigned.copy()

    def _getHashValue(self, attr):
        if attr in self._attributes:
            return unicode(self._attributes[attr] or u"")
        if attr in self._user:
            return unicode(self._user[attr] or u"")
        return u""


    def _fireSubscriptions(self, experiment, result):
        prev = self._assigned.get(experiment.key, None)
        if (
            not prev
            or prev[u"result"].inExperiment != result.inExperiment
            or prev[u"result"].variationId != result.variationId
        ):
            self._assigned[experiment.key] = {
                u"experiment": experiment,
                u"result": result,
            }
            for cb in self._subscriptions:
                try:
                    cb(experiment, result)
                except Exception:
                    pass


    def run(self, experiment):
        result = self._run(experiment)
        self._fireSubscriptions(experiment, result)
        return result

    def subscribe(self, callback):
        self._subscriptions.add(callback)
        return lambda: self._subscriptions.remove(callback)

    def _run(self, experiment, featureId = None):
        # 1. If experiment has less than 2 variations, return immediately
        if len(experiment.variations) < 2:
            return self._getExperimentResult(experiment, featureId=featureId)
        # 2. If growthbook is disabled, return immediately
        if not self._enabled:
            return self._getExperimentResult(experiment, featureId=featureId)
        # 2.5. If the experiment props have been overridden, merge them in
        if self._overrides.get(experiment.key, None):
            experiment.update(self._overrides[experiment.key])
        # 3. If experiment is forced via a querystring in the url
        qs = getQueryStringOverride(
            experiment.key, self._url, len(experiment.variations)
        )
        if qs is not None:
            return self._getExperimentResult(experiment, qs, featureId=featureId)
        # 4. If variation is forced in the context
        if self._forcedVariations.get(experiment.key, None) is not None:
            return self._getExperimentResult(
                experiment, self._forcedVariations[experiment.key], featureId=featureId
            )
        # 5. If experiment is a draft or not active, return immediately
        if experiment.status == u"draft" or not experiment.active:
            return self._getExperimentResult(experiment, featureId=featureId)
        # 6. Get the user hash attribute and value
        hashAttribute = experiment.hashAttribute or u"id"
        hashValue = self._getHashValue(hashAttribute)
        if not hashValue:
            return self._getExperimentResult(experiment, featureId=featureId)

        # 7. Exclude if user not in experiment.namespace
        if experiment.namespace and not inNamespace(hashValue, experiment.namespace):
            return self._getExperimentResult(experiment, featureId=featureId)

        # 7.5. If experiment has an include property
        if experiment.include:
            try:
                if not experiment.include():
                    return self._getExperimentResult(experiment, featureId=featureId)
            except Exception:
                return self._getExperimentResult(experiment, featureId=featureId)

        # 8. Exclude if condition is false
        if experiment.condition and not evalCondition(
            self._attributes, experiment.condition
        ):
            return self._getExperimentResult(experiment, featureId=featureId)

        # 8.1. Make sure user is in a matching group
        if experiment.groups and len(experiment.groups):
            expGroups = self._groups or {}
            matched = False
            for group in experiment.groups:
                if expGroups[group]:
                    matched = True
            if not matched:
                return self._getExperimentResult(experiment, featureId=featureId)
        # 8.2. If experiment.url is set, see if it's valid
        if experiment.url:
            if not self._urlIsValid(experiment.url):
                return self._getExperimentResult(experiment, featureId=featureId)

        # 9. Get bucket ranges and choose variation
        ranges = getBucketRanges(
            len(experiment.variations), experiment.coverage or 1, experiment.weights
        )
        n = gbhash(hashValue + experiment.key)
        assigned = chooseVariation(n, ranges)

        # 10. Return if not in experiment
        if assigned < 0:
            return self._getExperimentResult(experiment, featureId=featureId)

        # 11. If experiment is forced, return immediately
        if experiment.force is not None:
            return self._getExperimentResult(experiment, experiment.force, featureId=featureId)

        # 12. Exclude if in QA mode
        if self._qaMode:
            return self._getExperimentResult(experiment, featureId=featureId)

        # 12.5. If experiment is stopped, return immediately
        if experiment.status == u"stopped":
            return self._getExperimentResult(experiment, featureId=featureId)

        # 13. Build the result object
        result = self._getExperimentResult(experiment, assigned, True, featureId=featureId)

        # 14. Fire the tracking callback if set
        self._track(experiment, result)

        # 15. Return the result
        return result

    def _track(self, experiment, result):
        if not self._trackingCallback:
            return None
        key = (
            result.hashAttribute
            + unicode(result.hashValue)
            + experiment.key
            + unicode(result.variationId)
        )
        if not self._tracked.get(key):
            try:
                self._trackingCallback(experiment=experiment, result=result)
                self._tracked[key] = True
            except Exception:
                pass

    def _urlIsValid(self, pattern):
        if not self._url:
            return False

        try:
            r = re.compile(pattern)
            if r.search(self._url):
                return True

            pathOnly = re.sub(ur"^[^/]*/", u"/", re.sub(ur"^https?:\/\/", u"", self._url))
            if r.search(pathOnly):
                return True
            return False
        except Exception:
            return True

    def _getExperimentResult(
        self, experiment, variationId = -1, hashUsed = False, featureId = None
    ):
        hashAttribute = experiment.hashAttribute or u"id"

        inExperiment = True
        if variationId < 0 or variationId > len(experiment.variations) - 1:
            variationId = 0
            inExperiment = False

        return Result(
            featureId=featureId,
            inExperiment=inExperiment,
            variationId=variationId,
            value=experiment.variations[variationId],
            hashUsed=hashUsed,
            hashAttribute=hashAttribute,
            hashValue=self._getHashValue(hashAttribute),
        )
